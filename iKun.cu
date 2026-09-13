/**
 * Problem:
 *   There are KUN_NUMBER = 1000000 chickens (iKuns). Each chicken initially
 *   has two legs. In each operation, we uniformly pick one chicken that still
 *   has at least one leg, and cut off one of its legs. We repeat this
 *   operation KUN_NUMBER = 1000000 times. The question is: under mathematical
 *   expectation, how many chickens still have both legs intact at the end?
 *
 *   The theoretical expected value is approximately 317844.285.
 *
 * Simulation:
 *   Each GPU thread runs an independent Monte Carlo experiment. The result of
 *   each experiment is the number of intact two-legged chickens remaining
 *   after KUN_NUMBER operations. Thread results are summed on the device with
 *   atomicAdd and divided by the total number of experiments to estimate the
 *   expectation.
 *
 * Build and run:
 *   nvcc -O3 -o iKunCUDA iKun.cu
 *   ./iKunCUDA
 * 
 * Simulation Result: 
 *   E=317844.2806598758(158895695685787/499916800)
 *   Calculate Time: 6709569 milliseconds
 *
 *
 * Notes:
 *   - PCG32.h is CUDA-compatible. All PCG32 functions used here are decorated
 *     with PCG32_HOST_DEVICE and can be called inside CUDA kernels.
 *   - Each GPU thread must own its own PCG32Struct to avoid data races and to
 *     keep the generator state in registers/local memory.
 *   - The kernel uses PCG32Uniform_MaxBiggerThanMin(), which assumes min <= max
 *     and is faster than the general PCG32Uniform().
 */
#include <stdbool.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include "PCG32.h"
#include <vector>

/* CUDA launch configuration: 200 blocks * 512 threads = 102400 experiments. */
#define BLOCK_SIZE 512
#define GRID__SIZE 200

/* Number of chickens and number of operations per Monte Carlo experiment. */
#define KUN_NUMBER 1000000

/*
 * RepeatTimes is stored in constant memory so all GPU threads can read it
 * efficiently. RepeatTimesHost is the host-side copy used for the final
 * normalization.
 *
 * RepeatTimes is the total number of Monte Carlo experiments to run across
 * all threads, split evenly by the kernel.
 */
__constant__ long long unsigned int RepeatTimes=500000000LLU;
const long long unsigned int RepeatTimesHost=   500000000LLU;

/*
 * Device-side accumulator for the sum of all per-thread results.
 * Updated with atomicAdd to avoid write conflicts between threads.
 */
__device__ long long unsigned int ThreadResultDevice=0;

/*
 * CUDA kernel: each GPU thread runs one or more independent Monte Carlo
 * experiments.
 *
 * Parameters:
 *   status - device array of PCG32Struct, one generator per thread.
 *
 * Steps for each experiment:
 *   1. Start with KUN_NUMBER chickens that still have two legs
 *      (Kun2Legs = KUN_NUMBER), and KUN_NUMBER-1 chickens that still have
 *      one leg available to be cut (Kun1or2Leg = KUN_NUMBER-1).
 *   2. Perform KUN_NUMBER operations. In each operation, pick a random
 *      integer in [0, Kun1or2Leg]. If the value is less than Kun2Legs, it
 *      means we picked a two-legged chicken, so Kun2Legs decreases by 1.
 *      Otherwise we picked a one-legged chicken, so Kun1or2Leg decreases by 1.
 *   3. After all operations, Kun2Legs is the number of chickens that still
 *      have both legs. Add it to myResult.
 *
 * After all experiments assigned to this thread, the accumulated result is
 * added to the global device counter with atomicAdd.
 *
 * PCG32 API used:
 *   PCG32Uniform_MaxBiggerThanMin(PCG32Struct* status,
 *                                 unsigned min,
 *                                 unsigned max)
 *     Returns a uniform integer in [min, max]. It assumes min <= max and is
 *     faster than PCG32Uniform() because it skips the min/max ordering check.
 *
 *   PCG32SetMultipleSeeds() is called on the host to initialize all generators
 *   before the kernel is launched.
 */
__global__ void threadSimulate(PCG32Struct* status){
    long long unsigned int myThreadIndex=blockIdx.x*blockDim.x+threadIdx.x;

    /*
     * Copy this thread's generator state from global memory to local memory.
     * Local memory / registers are much faster than repeated global accesses.
     */
    PCG32Struct PCGStatus=status[myThreadIndex];

    /* Number of Monte Carlo experiments assigned to this thread. */
    const long long unsigned int myRepeatTimes=(RepeatTimes/BLOCK_SIZE)/GRID__SIZE;
    long long unsigned int myResult=0;

    for(unsigned indexRepeat=0;indexRepeat<myRepeatTimes;indexRepeat=indexRepeat+1){
        /*
         * State of one Monte Carlo experiment:
         *   Kun2Legs   - number of chickens that still have two legs.
         *   Kun1or2Leg - number of chickens that still have at least one leg
         *                available to be cut, minus 1. This is used as the
         *                inclusive upper bound for random selection.
         */
        unsigned Kun2Legs=KUN_NUMBER;
        unsigned Kun1or2Leg=KUN_NUMBER-1;

        /*
         * Perform KUN_NUMBER operations. Each operation picks one chicken
         * uniformly from those that still have at least one leg.
         */
        for(unsigned indexOperation=0;indexOperation<KUN_NUMBER;indexOperation=indexOperation+1){
            if(PCG32Uniform_MaxBiggerThanMin(&PCGStatus,0,Kun1or2Leg)<Kun2Legs){
                /* Picked a two-legged chicken: it loses one leg. */
                Kun2Legs=Kun2Legs-1;
            }else{
                /* Picked a one-legged chicken: it loses its last leg. */
                Kun1or2Leg=Kun1or2Leg-1;
            }
        }

        /* Kun2Legs is the number of intact two-legged chickens in this run. */
        myResult=myResult+Kun2Legs;
    }

    /*
     * Accumulate this thread's total into the global device counter.
     * atomicAdd ensures correctness when many threads finish simultaneously.
     */
    atomicAdd(&ThreadResultDevice,myResult);
}

int main(int argc, char const *argv[]){
	struct timeval start,end;
    unsigned milliseconds=0;
    gettimeofday(&start,NULL);

    /*
     * Generate distinct seeds for all GPU threads. The base seed comes from
     * the current time, and each thread gets a different offset.
     */
    long long unsigned int BaseSeed=PCG32TimeNanoeconds();
    std::vector<PCG32Struct> hostStatus(GRID__SIZE*BLOCK_SIZE);
    std::vector<long long unsigned int> hostSeeds(hostStatus.size());
    for(unsigned threadIndex=0;threadIndex<hostStatus.size();threadIndex=threadIndex+1){
        hostSeeds[threadIndex]=BaseSeed+threadIndex;
    }

    /*
     * Initialize all PCG32 generators on the host.
     *
     * PCG32SetMultipleSeeds(PCG32Struct* statusArray,
     *                       const long long unsigned int* baseSeeds,
     *                       unsigned count)
     *
     * Each generator gets its own seed and is advanced by a sequence of
     * primes to decorrelate the random streams.
     */
    PCG32SetMultipleSeeds(hostStatus.data(),hostSeeds.data(),hostStatus.size());

    /*
     * Copy the initialized generators to the device.
     * Each GPU thread will read its own PCG32Struct from this array.
     */
    PCG32Struct* deviceStatus;
    cudaMalloc((void**)(&deviceStatus),hostStatus.size()*sizeof(PCG32Struct));
    cudaMemcpy(deviceStatus,hostStatus.data(),hostStatus.size()*sizeof(PCG32Struct),cudaMemcpyHostToDevice);

    /* Launch the Monte Carlo kernel. */
    threadSimulate<<<GRID__SIZE,BLOCK_SIZE>>>(deviceStatus);

    cudaDeviceSynchronize();

    /*
     * Retrieve the global sum of all experiments from device memory.
     * ThreadResultDevice lives in device global memory, so we use
     * cudaMemcpyFromSymbol.
     */
    long long unsigned int ThreadResult;
    cudaMemcpyFromSymbol(&ThreadResult,ThreadResultDevice,sizeof(long long unsigned int));
    cudaFree(deviceStatus);

    /*
     * Compute the actual number of Monte Carlo experiments that were run.
     * The kernel divides RepeatTimes by BLOCK_SIZE and GRID__SIZE using
     * integer division, so the effective count may be slightly smaller than
     * RepeatTimesHost.
     */
    long long unsigned int RealRepeatTimes=((RepeatTimesHost/BLOCK_SIZE)/GRID__SIZE)*BLOCK_SIZE*GRID__SIZE;
    printf("E=%.10lf(%llu/%llu)\n",(double)ThreadResult/(double)RealRepeatTimes,ThreadResult,RealRepeatTimes);

    gettimeofday(&end,NULL);
    milliseconds=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0+0.5;
    printf("Calculate Time: %u milliseconds\n",milliseconds);
    
	return 0;
}