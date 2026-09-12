/**
 * This program estimates the average projected area of a randomly oriented
 * unit cube when sunlight is perpendicular to the ground.
 *
 * The theoretical expected projected area of a unit cube under uniform random
 * rotation is exactly 1.5. This program uses Monte Carlo simulation to verify
 * that value.
 *
 * Two implementations are compared:
 *   1. Scalar version:
 *      Each OpenMP thread owns one PCG32Struct and generates one standard
 *      normal sample at a time using PCG32StandardNormal().
 *   2. AVX512 version:
 *      Each OpenMP thread owns one __x16__PCG32Struct and generates 16
 *      standard normal samples at a time using __x16__PCG32StandardNormal().
 *
 * Build:
 *   clear && gcc -o CubeProjection CubeProjection.c -O2 -mavx512f -mavx512dq -mfma -lm -fopenmp && ./CubeProjection
 *
 * Notes:
 *   - AVX512 path requires GCC on x86_64 with AVX512F, AVX512DQ and FMA.
 *   - OpenMP is used to parallelize Monte Carlo sampling across threads.
 *   - Each thread must use its own PCG32 generator to avoid data races.
 * 
 * Simulation Result: 
 * 
 * AVX512 version
 * average projected area: 1.4999637668
 * average projected area: 1.5000020164
 * average projected area: 1.4999914889
 * average projected area: 1.4999829889
 * average projected area: 1.5000011244
 * Calculate Prime Time: 137 milliseconds
 * scalar version
 * average projected area: 1.5000223160
 * average projected area: 1.5000080209
 * average projected area: 1.5000065327
 * average projected area: 1.4999800102
 * average projected area: 1.4999847645
 * Calculate Prime Time: 608 milliseconds
 */
#include <stdbool.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include "PCG32.h"  // #include "PCG32.h"
#include <string.h>
#include <omp.h>
#include <math.h>
#include <time.h>

/* Number of OpenMP threads used by both scalar and AVX512 simulations. */
#define THREAD_NUMBER 32

/*
 * Total number of Monte Carlo samples per timing run.
 * Scalar: each thread handles RepeatTimes / THREAD_NUMBER samples.
 * AVX512: each batch produces 16 samples, so each thread 
 * handles RepeatTimes / THREAD_NUMBER / 16 batches.
 */
const long long unsigned int RepeatTimes=100000000LLU;

/* Pi used for geometric calculations. */
const double Pi=3.1415926535897932384626433832795028841971693993751;

/*
 * Compute the Euclidean norm of a 4D vector.
 *
 * A random rotation is represented as a quaternion. Four independent standard
 * normal variables form a 4D Gaussian vector, which is normalized into a
 * uniformly distributed unit quaternion on S^3.
 */
static inline double Normalize(const double vector[4]){
    return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]+vector[3]*vector[3]);
}

/*
 * Compute the projected area of a unit cube for one random quaternion.
 *
 * Input random[4] contains four independent standard normal variables.
 * They are normalized into a unit quaternion q=(q0,q1,q2,q3).
 *
 * The third column of the rotation matrix represented by q is the rotated
 * z-axis, i.e. the normal vector of one cube face:
 *
 *   n_x = 2*(q1*q3-q0*q2)
 *   n_y = 2*(q2*q3+q0*q1)
 *   n_z = 1-2*(q1*q1+q2*q2)
 *
 * For a unit cube, the projected area onto the ground plane is the L1 norm
 * of this unit normal:
 *
 *   area = |n_x|+|n_y|+|n_z|
 *
 * Its expectation over uniform random rotations is exactly 1.5.
 */
static inline double Projection(const double random[4]){
    const double distance=Normalize(random);
    const double q[4]={
        random[0]/distance,
        random[1]/distance,
        random[2]/distance,
        random[3]/distance
    };
    const double Normal[3]={
        (q[1]*q[3]-q[0]*q[2])*2,
        (q[2]*q[3]+q[0]*q[1])*2,
        1.0-(q[1]*q[1]+q[2]*q[2])*2
    };
    return fabs(Normal[0])+fabs(Normal[1])+fabs(Normal[2]);
}

/*
 * Scalar Monte Carlo simulation for one OpenMP thread.
 *
 * Each thread owns one PCG32Struct in status[]. The thread generates four
 * standard normal samples per Monte Carlo sample using
 * PCG32StandardNormal(PCG32Struct* status).
 *
 * PCG32StandardNormal returns one N(0,1) sample by the Box-Muller method.
 * The generator state is kept in the caller-owned PCG32Struct, so each
 * thread must use its own structure.
 *
 * The average projected area for this thread is written to
 * threadResult[omp_get_thread_num()].
 */
void ThreadSimulateScalar(PCG32Struct* status,double* threadResult){
    const long long unsigned int myThreadIndex=omp_get_thread_num();
    double myResult=0.0;

    /* Divide total samples evenly across OpenMP threads. */
    const long long unsigned int myRepeatTimes=RepeatTimes/omp_get_num_threads();

    /* Each thread uses its own scalar generator to avoid data races. */
    PCG32Struct* myStatus=status+myThreadIndex;
    for(long long unsigned int times=0;times<myRepeatTimes;times=times+1){
        /* Four independent standard normal variables form one quaternion. */
        const double random[4]={PCG32StandardNormal(myStatus),PCG32StandardNormal(myStatus),PCG32StandardNormal(myStatus),PCG32StandardNormal(myStatus)};
        myResult=myResult+Projection(random);
    }

    /* Store this thread's average projected area. */
    myResult=myResult/myRepeatTimes;
    threadResult[myThreadIndex]=myResult;
}

/*
 * AVX512 Monte Carlo simulation for one OpenMP thread.
 *
 * Each thread owns one __x16__PCG32Struct. The AVX512 API generates 16
 * standard normal samples per call. Four calls produce 4*16=64 doubles,
 * which are interpreted as 16 quaternions of 4 components each.
 *
 * AVX512 PCG32 API used:
 *   __x16__PCG32StandardNormal(__x16__PCG32Struct* status,
 *                              __x16__DoubleArray random)
 *     Fills random[0..15] with 16 independent standard normal samples.
 *
 * Data layout:
 *   __x16__DoubleArray randoms[4] is a 4x16 array of doubles.
 *   Casting it to double* gives 64 contiguous doubles:
 *     randomArray[0..3]   -> quaternion 0
 *     randomArray[4..7]   -> quaternion 1
 *     ...
 *     randomArray[60..63] -> quaternion 15
 */
void ThreadSimulateAVX512(__x16__PCG32Struct* AVX512Status,double* threadResult){
    const long long unsigned int myThreadIndex=omp_get_thread_num();
    double myResult=0.0;

    /* Four AVX512 batches, each holding 16 doubles. */
    __x16__DoubleArray randoms[4];

    /* Flat view: 4*16=64 contiguous doubles, i.e. 16 quaternions. */
    double* randomArray=(double*)randoms;

    /* Each batch produces 16 samples. */
    const long long unsigned int myRepeatTimes=RepeatTimes/omp_get_num_threads()/16;

    /* Each thread uses its own AVX512 generator. */
    __x16__PCG32Struct* myStatus=AVX512Status+myThreadIndex;
    for(long long unsigned int times=0;times<myRepeatTimes;times=times+1){
        /*
         * Generate 4*16 standard normal variables.
         * Every 4 consecutive doubles form one quaternion.
         */
        __x16__PCG32StandardNormal(myStatus,randoms[0]);
        __x16__PCG32StandardNormal(myStatus,randoms[1]);
        __x16__PCG32StandardNormal(myStatus,randoms[2]);
        __x16__PCG32StandardNormal(myStatus,randoms[3]);
        double __x16__result=0.0;
        for(unsigned index=0;index<4*16;index=index+4){
            myResult=myResult+Projection(randomArray+index);
        }
    }

    /* Normalize by number of batches and by 16 samples per batch. */
    myResult=myResult/myRepeatTimes/16;
    threadResult[myThreadIndex]=myResult;
}

int main(int argc,char const *argv[]){
    struct timeval start,end;
    unsigned milliseconds=0;

    /* Per-thread simulation results. */
    double result[THREAD_NUMBER];
    double averageResult;

    /* Number of repeated timing runs. */
    const unsigned repeat=5;

    /*
     * Seeds for scalar generators: one 64-bit seed per thread.
     * Seeds for AVX512 generators: one __x16__SeedArray per thread,
     * where each __x16__SeedArray contains 16 64-bit lane seeds.
     */
    long long unsigned int seeds[THREAD_NUMBER];
    __x16__SeedArray AVX512Seeds[THREAD_NUMBER];

    /*
     * Generator states:
     *   status[]       - one scalar PCG32Struct per thread.
     *   AVX512Status[] - one AVX512 PCG32Struct per thread, 16 lanes each.
     */
    PCG32Struct status[THREAD_NUMBER];
    __x16__PCG32Struct AVX512Status[THREAD_NUMBER];

    /*
     * Build time-based seeds with per-thread offsets to reduce correlation
     * between threads and SIMD lanes.
     */
    long long unsigned int time=PCG32TimeNanoeconds();
    for(unsigned threadIndex=0;threadIndex<THREAD_NUMBER;threadIndex=threadIndex+1){
        seeds[threadIndex]=time+threadIndex*0x123456789ABCDEFLLU+1123;
        for(unsigned index=0;index<16;index=index+1){
            AVX512Seeds[threadIndex][index]=time+threadIndex*0x123456789ABCDEFLLU+index;
        }
    }

    /*
     * Initialize all scalar generators at once.
     *
     * PCG32SetMultipleSeeds(PCG32Struct* statusArray,
     *                       const long long unsigned int* baseSeeds,
     *                       unsigned count)
     *
     * Each generator gets its own base seed and is advanced by a sequence
     * of primes to further decorrelate the streams.
     */
    PCG32SetMultipleSeeds(status,seeds,THREAD_NUMBER);

    /*
     * Initialize all AVX512 generators at once.
     *
     * __x16__PCG32SetMultipleSeeds(__x16__PCG32Struct* status,
     *                              const __x16__SeedArray baseSeed[],
     *                              unsigned count)
     *
     * Each AVX512 generator contains 16 independent PCG32 lanes.
     * All lanes of all generators are initialized from AVX512Seeds.
     */
    __x16__PCG32SetMultipleSeeds(AVX512Status,AVX512Seeds,THREAD_NUMBER);

    /* ---------------- AVX512 version ---------------- */
    printf("AVX512 version\n");
    gettimeofday(&start,NULL);
    for(unsigned index=0;index<repeat;index=index+1){
        /*
         * OpenMP parallel region. Each thread calls ThreadSimulateAVX512()
         * and uses its own AVX512 generator.
         */
        #pragma omp parallel
        {
            ThreadSimulateAVX512(AVX512Status,result);
        }
        averageResult=0.0;
        for(unsigned index=0;index<THREAD_NUMBER;index=index+1){
            averageResult=averageResult+result[index];
        }
        averageResult=averageResult/THREAD_NUMBER;
        printf("average projected area: %.10lf\n",averageResult);
    }
    gettimeofday(&end,NULL);
    milliseconds=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0+0.5;
    printf("Calculate Prime Time: %u milliseconds\n",milliseconds);

    /* ---------------- scalar version ---------------- */
    printf("scalar version\n");
    gettimeofday(&start,NULL);
    for(unsigned index=0;index<repeat;index=index+1){
        /*
         * OpenMP parallel region. Each thread calls ThreadSimulateScalar()
         * and uses its own scalar generator.
         */
        #pragma omp parallel
        {
            ThreadSimulateScalar(status,result);
        }
        averageResult=0.0;
        for(unsigned index=0;index<THREAD_NUMBER;index=index+1){
            averageResult=averageResult+result[index];
        }
        averageResult=averageResult/THREAD_NUMBER;
        printf("average projected area: %.10lf\n",averageResult);
    }
    gettimeofday(&end,NULL);
    milliseconds=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0+0.5;
    printf("Calculate Prime Time: %u milliseconds\n",milliseconds);

    return 0;
}