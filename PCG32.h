/**
// 
// PCG32.h
// High-performance PCG-XSH-RR pseudorandom number generator with comprehensive distributions
//
// Lightweight, high-performance pseudorandom number generator based on the PCG-XSH-RR algorithm,
// offering excellent statistical quality and computational efficiency. Provides both C and C++
// interfaces with full CUDA compatibility for GPU acceleration. Some function provides AVX512 
// version. 
// 
// @author Pandora Eartha
// @see https://github.com/PandoraEartha/HighPerformanceRandomEngine
//
// 
// C API (prefix: PCG32)
//
// | Function                                    | Description                                   |
// |---------------------------------------------|-----------------------------------------------|
// | PCG32SetSingleSeed                          | Initialize generator with single seed         |
// | PCG32SetMultipleSeeds                       | Initialize multiple generators with seeds     |
// | PCG32                                       | Generate 32-bit uniform integer               |
// | PCG32Uniform                                | Uniform integer [min, max] (auto-order)       |
// | PCG32Uniform_Strict                         | Uniform integer (min <= max, gap not power-2) |
// | PCG32UniformSetStrictRange                  | Preset strict range for repeated use          |
// | PCG32Uniform_StrictRangeUnchanged           | Generate using preset strict range            |
// | PCG32Uniform_MaxBiggerThanMin               | Fast uniform (assumes min <= max)             |
// | PCG32UniformReal                            | Uniform real in [min, max)                    |
// | PCG32StandardNormal                         | Standard normal N(0,1)                        |
// | PCG32GammaInitialize                        | Initialize Gamma distribution (alpha >= 1)    |
// | PCG32Gamma                                  | Generate Gamma(alpha, beta)                   |
// | PCG32Binomial                               | Binomial(n, p)                                |
// | PCG32PoissonInitialize                      | Initialize Poisson distribution               |
// | PCG32Poisson                                | Generate Poisson(mu)                          |
// | PCG32Exponential                            | Exponential(lambda)                           |
// | PCG32PowerLaw                               | Power-law (min, alpha)                        |
// | PCG32Geometric                              | Geometric(p)                                  |
// | PCG32Geometric_SmallProbability             | Geometric(p) optimized for small p            |
// | PCG32LogNormal                              | Log-normal(mu, sigma)                         |
// | PCG32StandardLogNormal                      | Standard log-normal (mu=0, sigma=1)           |
// | PCG32Benford                                | Benford's law (digits 1-12)                   |
// | PCG32Benford_SpecificLength                 | Benford with specified digit length range     |
// | PCG32RandomPointInSphere3D                  | Uniform point in 3D sphere                    |
// | PCG32RandomPointInSphereNDimension          | Uniform point in N-dimensional sphere         |
// | PCG32RandomPointInCycle                     | Uniform point in 2D circle                    |
// | PCG32UniformSumReal                         | N uniform reals summing to a fixed value      |
// | PCG32UniformShuffle                         | Fisher-Yates shuffle (macro, type-generic)    |
// | PCG32UniformShuffle_FirstK                  | Fisher-Yates shuffle first K element          |
// | PCG32TimeNanoeconds                         | Get current time in nanoseconds since epoch   |
// 
// Features:
// - 32-bit uniform random integers [0, 0xFFFFFFFF]
// - Uniform real numbers [min, max)
// - Multiple integer uniform variants (general, strict, fast, preset-range)
// - Standard normal distribution (Box-Muller transform)
// - Gamma distribution (alpha >= 1, Marsaglia & Tsang method)
// - Binomial distribution (efficient BTPE algorithm)
// - Poisson distribution (saddlepoint approximation for large mu)
// - Exponential distribution
// - Power-law (Pareto) distribution
// - Geometric distribution (standard + small-probability optimized)
// - Log-normal distribution (standard + general)
// - Benford's law distributed numbers
// - N-dimensional uniform points in spheres
// - Uniform points on/inside 2D circles
// - Uniform simplex sampling (fixed-sum real variables)
// - Fisher-Yates array shuffling (template for arbitrary types)
// - Time-based seed generation via PCG32TimeNanoeconds
//
// AVX512 Vectorization Support
//
// When compiled with GCC or G++ on x86_64 platforms with AVX512F and AVX512DQ instruction sets,
// the library provides SIMD functions that generate 16 random numbers simultaneously.
//
// Build flags: -mavx512f -mavx512dq -mfma
//
// AVX512 data types (64-byte aligned):
//   __x16__StateArray    : 16 x uint64_t (generator states)
//   __x16__SeedArray     : 16 x uint64_t (seeds)
//   __x16__UnsignedArray : 16 x uint32_t (random integers)
//   __x16__DoubleArray   : 16 x double   (random reals)
//
// AVX512 functions:
//   __x16__PCG32SetSingleSeed                 : Initialize 16 generators with a seed array
//   __x16__PCG32SetMultipleSeeds              : Initialize multiple 16-lane generator sets
//   __x16__PCG32                              : Generate 16 random 32-bit integers
//   __x16__PCG32UniformReal                   : Generate 16 uniform reals in [0,1)
//   __x16__PCG32UniformReal_MinMax            : Generate 16 uniform reals in [min, max)
//   __x16__PCG32UniformSetStrictRange         : Preset a strict range for all 16 lanes
//   __x16__PCG32Uniform_StrictRangeUnchanged  : Generate 16 integers from preset range
//   __x16__PCG32StandardNormal                : Generate 16 standard normal N(0,1) samples
//   __x16__PCG32Normal                        : Generate 16 normal N(mu, sigma) samples
//
//
// C++ Class API (PCG32PRNG)
//
// All C functions are wrapped as methods of the PCG32PRNG class. Template support
// is provided for UniformShuffle and UniformShuffle_FirstK.
//
// CUDA Support
//
// All functions are decorated with PCG32_HOST_DEVICE, enabling direct usage within
// CUDA kernel code. Simply include this header in your .cu files.
//
// Important Notes
//
// 1. Seed required: Always call PCG32SetSingleSeed/PCG32SetMultipleSeeds before generating numbers.
// 2. Multi-threading: Use different PCG32Struct per thread with different seeds.
// 3. Gamma restriction: Currently only supports shape parameter alpha >= 1.
// 4. Performance tip: For repeated generation within the same integer range,
//    use UniformSetStrictRange + Uniform_StrictRangeUnchanged for optimal speed.
// 5. AVX512 restriction: AVX512 functions are only available with GCC/G++ on x86_64
//    and require -mavx512f -mavx512dq  -mfma flags. They are not available in CUDA mode.
//
// Usage Examples
//
// C++ Interface:
//
//   #include "PCG32.h"
//   
//   // Basic initialization and uniform real
//   PCG32PRNG rng(PCG32TimeNanoeconds(nullptr));
//   double u = rng.UniformReal(-1.0, 1.0);  // uniform in [-1.0, 1.0)
//   
//   // Binomial distribution: 1000 trials with p=0.9
//   unsigned binom = rng.Binomial(0.9, 1000);
//   
//   // Gamma distribution: must initialize first (alpha >= 1)
//   rng.GammaInitialize(2.0, 1.0))           // alpha=2.0, beta=1.0
//   double gamma = rng.Gamma();              // Gamma(2.0, 1.0) sample
//   
//   // Fast strict-range integer generation: preset range once, generate many times
//   rng.UniformSetStrictRange(0, 999);       // preset range [0, 999]
//   for (int i = 0; i < 100; ++i) {
//       unsigned x = rng.Uniform_StrictRangeUnchanged();  // fast, no range checks
//   }
//   
//   // Random integer with auto-range (min/max auto-sorted)
//   unsigned x = rng.Uniform(10, 20);
//   
//   // Shuffle an array
//   int myArray[100];
//   rng.UniformShuffle(myArray, 100);
//
// C Interface:
//
//   #include "PCG32.h"
//   
//   PCG32Struct state;
//   PCG32SetSingleSeed(&state, PCG32TimeNanoeconds(NULL));
//   
//   // Uniform real in [-1.0, 1.0)
//   double u = PCG32UniformReal(&state, -1.0, 1.0);
//   
//   // Binomial: 1000 trials with p=0.9
//   unsigned binom = PCG32Binomial(&state, 0.9, 1000);
//   
//   // Gamma: must initialize first
//   if (PCG32GammaInitialize(&state, 2.0, 1.0)) {
//       double gamma = PCG32Gamma(&state);
//   }
//   
//   // Fast strict-range preset
//   PCG32UniformSetStrictRange(&state, 0, 999);
//   for (int i = 0; i < 100; ++i) {
//       unsigned x = PCG32Uniform_StrictRangeUnchanged(&state);
//   }
//   
//   // Basic random integer
//   unsigned x = PCG32Uniform(&state, 10, 20);
//   
//   // Shuffle array (macro, type-generic)
//   int myArray[100];
//   PCG32UniformShuffle(&state, myArray, 100);
//
// Multi-threaded Seed Initialization with PCG32SetMultipleSeeds:
//
//   #include "PCG32.h"
//   #include <omp.h>
//   
//   #define THREAD_NUMBER 32
//
//   // or use malloc or std::vector if THREAD_NUMBER is not a compile-time constant
//   long long unsigned int seeds[THREAD_NUMBER]; 
//   PCG32Struct status[THREAD_NUMBER];
//   long long unsigned int time = PCG32TimeNanoeconds();
//   
//   // Generate distinct seeds for each thread
//   for (unsigned threadIndex = 0; threadIndex < THREAD_NUMBER; threadIndex++) {
//       seeds[threadIndex] = time + threadIndex * 0x123456789ABCDEFLLU + 1123;
//   }
//   
//   // Initialize all generators with one call
//   PCG32SetMultipleSeeds(status, seeds, THREAD_NUMBER);
//   
//   #pragma omp parallel num_threads(THREAD_NUMBER)
//   {
//       int tid = omp_get_thread_num();
//       // Each thread uses its own generator
//       double u = PCG32UniformReal(&status[tid], 0.0, 1.0);
//       // ... generate more numbers ...
//   }
//
// CUDA Seed Initialization with PCG32SetMultipleSeeds:

//   #include "PCG32.h"
//   
//   // CUDA kernel - each thread processes its own generator
//   __global__ void kernel(PCG32Struct* deviceStatus, unsigned* results, int N) {
//       int idx = blockIdx.x * blockDim.x + threadIdx.x;
//       if (idx >= N) return;
//       
//       // Each thread copies its status from device memory to local stack
//       // (register/local memory is much faster than global memory)
//       PCG32Struct localStatus = deviceStatus[idx];
//       
//       // Generate random numbers using local state
//       for (int i = 0; i < 10; i++) {
//           unsigned random = PCG32(&localStatus);
//           results[idx * 10 + i] = random;
//       }
//       
//       // Write back updated state if needed for subsequent kernel launches
//       deviceStatus[idx] = localStatus;
//   }
//   
//   int main() {
//       const int N = 1024 * 1024;  // 1 million threads
//       const int BLOCK_SIZE = 256;
//       const int GRID_SIZE = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
//       
//       // Host-side arrays
//       long long unsigned int* seeds = (long long unsigned int*)malloc(N * sizeof(long long unsigned int));
//       PCG32Struct* hostStatus = (PCG32Struct*)malloc(N * sizeof(PCG32Struct));
//       unsigned* hostResults = (unsigned*)malloc(N * 10 * sizeof(unsigned));
//       
//       // Device-side arrays
//       PCG32Struct* deviceStatus;
//       unsigned* deviceResults;
//       cudaMalloc(&deviceStatus, N * sizeof(PCG32Struct));
//       cudaMalloc(&deviceResults, N * 10 * sizeof(unsigned));
//       
//       // Generate seeds on host
//       long long unsigned int time = PCG32TimeNanoeconds();
//       for (int i = 0; i < N; i++) {
//           seeds[i] = time + i * 0x123456789ABCDEFLLU + 1123;
//       }
//       
//       // Initialize all generators on host using PCG32SetMultipleSeeds
//       PCG32SetMultipleSeeds(hostStatus, seeds, N);
//       
//       // Copy initialized generators to device
//       cudaMemcpy(deviceStatus, hostStatus, N * sizeof(PCG32Struct), cudaMemcpyHostToDevice);
//       
//       // Launch kernel with proper grid configuration
//       kernel<<<GRID_SIZE, BLOCK_SIZE>>>(deviceStatus, deviceResults, N);
//       cudaDeviceSynchronize();
//       
//       // Copy results back
//       cudaMemcpy(hostResults, deviceResults, N * 10 * sizeof(unsigned), cudaMemcpyDeviceToHost);
//       
//       // Clean up
//       free(seeds);
//       free(hostStatus);
//       free(hostResults);
//       cudaFree(deviceStatus);
//       cudaFree(deviceResults);
//       
//       return 0;
//   }
//
// AVX512 Batch Generation (compile with -mavx512f -mavx512dq):
//
//   #include "PCG32.h"
//   
//   // 1. __x16__PCG32SetSingleSeed - Initialize 16 generators with seeds
//   __x16__PCG32Struct avxState;
//   __x16__SeedArray seeds = {
//       0x123456789ABCDEF0ULL, 0x23456789ABCDEF01ULL,  // ... 16 seeds total
//       // ... fill all 16 seeds
//   };
//   __x16__PCG32SetSingleSeed(&avxState, seeds);
//   
//   // 2. __x16__PCG32 - Generate 16 random 32-bit integers
//   __x16__UnsignedArray randoms;
//   __x16__PCG32(&avxState, randoms);
//   // randoms[0] through randoms[15] contain random 32-bit values
//   
//   // 3. __x16__PCG32UniformReal - Generate 16 uniform reals in [0, 1)
//   __x16__DoubleArray uniforms;
//   __x16__PCG32UniformReal(&avxState, uniforms);
//   // uniforms[0] through uniforms[15] are in [0.0, 1.0)
//   
//   // 4. __x16__PCG32UniformSetStrictRange - Preset strict range for all lanes
//   //    Range: [0, 999] for all 16 generators
//   __x16__PCG32UniformSetStrictRange(&avxState, 0, 999);
//   
//   // 5. __x16__PCG32Uniform_StrictRangeUnchanged - Generate 16 integers from preset range
//   __x16__UnsignedArray strictRandoms;
//   for (int batch = 0; batch < 10; ++batch) {
//       __x16__PCG32Uniform_StrictRangeUnchanged(&avxState, strictRandoms);
//       // strictRandoms[0..15] are all in [0, 999]
//       // Process the 16 random values...
//   }
//   // Note: For optimal performance with AVX512 strict range, the range should be
//   // preset once with __x16__PCG32UniformSetStrictRange, then reuse it with
//   // __x16__PCG32Uniform_StrictRangeUnchanged for many batches.
// 
//   // 6. Generate 16 uniform reals in custom range [min, max)
//   __x16__DoubleArray customUniforms;
//   __x16__PCG32UniformReal_MinMax(&avxState, -5.0, 5.0, customUniforms);
//   // customUniforms[0..15] are all in [-5.0, 5.0)
//
//   // 7. Generate 16 standard normal N(0,1) samples
//   __x16__DoubleArray normalSamples;
//   __x16__PCG32StandardNormal(&avxState, normalSamples);
//   // normalSamples[0..15] are standard normal distributed
//   
//   // 8. Generate 16 normal N(mu, sigma) samples
//   __x16__DoubleArray customNormalSamples;
//   __x16__PCG32Normal(&avxState, customNormalSamples, 2.0, 0.5);
//   // customNormalSamples[0..15] are N(2.0, 0.5) distributed
//
// AVX512 Multi-threaded Seed Initialization:
//
//   #include "PCG32.h"
//   #include <omp.h>
//
//   #define THREAD_NUMBER 32
//
//   // or use malloc or std::vector if THREAD_NUMBER is not a compile-time constant
//   __x16__SeedArray AVX512Seeds[THREAD_NUMBER];
//   __x16__PCG32Struct AVX512Status[THREAD_NUMBER];
//   long long unsigned int time = PCG32TimeNanoeconds();
//   
//   // Generate distinct seeds for each thread
//   for (unsigned threadIndex = 0; threadIndex < THREAD_NUMBER; threadIndex++) {
//       for (unsigned index = 0; index < 16; index++) {
//           AVX512Seeds[threadIndex][index] = time + threadIndex * 0x123456789ABCDEFLLU + index;
//       }
//   }
//   
//   // Initialize all generators for scalar and AVX512 usage
//   __x16__PCG32SetMultipleSeeds(AVX512Status, AVX512Seeds, THREAD_NUMBER);
//   
//   #pragma omp parallel num_threads(THREAD_NUMBER)
//   {
//       int tid = omp_get_thread_num();
//       
//       // AVX512 generator for this thread (generates 16 values at once)
//       __x16__DoubleArray avxUniforms;
//       __x16__PCG32UniformReal(&AVX512Status[tid], avxUniforms);
//       // Process 16 random values...
//   }
//
// Complete AVX512 compilation example:
//   g++ -O3 -mavx512f -mavx512dq -mfma -std=c++11 -o myapp main.cpp
// 
 */

#ifndef __PCG32_H__
#define __PCG32_H__

#if defined(__cplusplus)||PCG32_CUDA
    #define PCG32CXX 1
#else
    #define PCG32CXX 0
#endif

#if defined(__CUDACC__)||defined(__CUDA_ARCH__)||defined(__CUDA_LIBDEVICE__)
    #define PCG32_CUDA 1
#else
    #define PCG32_CUDA 0
#endif

#if defined(__GNUC__)&&defined(__x86_64__)&&defined(__AVX512F__)
    #define PCG32_AVX512 1
#else
    #define PCG32_AVX512 0
#endif

#if PCG32_CUDA
    #define PCG32_HOST_DEVICE __host__ __device__
    #define PCG32_DEVICE               __device__
    #define PCG32_HOST        __host__
    #include <cuda_runtime_api.h>
    #include <limits>
    #include <float.h>
    #include <time.h>
#else
    #define PCG32_HOST_DEVICE
    #define PCG32_DEVICE
    #define PCG32_HOST
    #include <math.h>
    #include <stdbool.h>
    #include <float.h>
    #include <stdio.h>
    #include <time.h>
#endif

#if PCG32_CUDA
    #define PCG32MULTIPLIER 6364136223846793005LLU
    #define PCG32INCREMENT  1442695040888963407LLU
    #define PCG32MAX        0x00000000FFFFFFFFLLU
    #define PCG32FALSE 0
    #define PCG32TRUE  1
    #define PCG32REAL_SCALE (1.0/(PCG32MAX+1))
#else
    static long long unsigned int const PCG32MULTIPLIER=6364136223846793005LLU;
    static long long unsigned int const PCG32INCREMENT =1442695040888963407LLU;
    static long long unsigned int const PCG32MAX       =0x00000000FFFFFFFFLLU;
    static const unsigned char PCG32FALSE=0;
    static const unsigned char PCG32TRUE =1;
    static const double PCG32REAL_SCALE=(double)1/(double)(PCG32MAX+1);
#endif

#if defined(__GNUC__)||defined(__clang__)
    #define PCG32_DEPRECATED(message) __attribute__((deprecated(message)))
#elif defined(_MSC_VER)
    #define PCG32_DEPRECATED(message) __declspec(deprecated(message))
#else
    #define PCG32_DEPRECATED(message)
#endif

#define PCG32_INT_MAX      0x7FFFFFFFU
#define PCG32_UNSIGNED_MAX 0xFFFFFFFFU

typedef struct PCG32Struct{
    long long unsigned int state;
    double normalDistributionSaved;
    double gammaD;
    double gammaC;
    double gammaBeta;
    double gammaNegativeSqrt9AlphaSub3;
    unsigned normalDistributionSavedValid;
    unsigned uniformStrictGap;
    unsigned uniformStrictRange;
    unsigned uniformStrictMin;
    struct PoissonParameter{
        double mu;
        double expmu;
        double floorMu;
        double dx;
        double roundMu;
        double logMu;
        double lgmma;
        double sqrtMu;
        double c[6];
        double cx;
        double sqrtCx;
        double cxDiviedBy1;
        double c2b;
        double cb;
        double cx2;
    }Poisson;
}PCG32Struct;


// Auxiliary C Functions
PCG32_HOST_DEVICE static inline unsigned rotr32(unsigned x,unsigned r);
PCG32_HOST_DEVICE static inline long long unsigned int PCG32NextPrime(long long unsigned int prime);
// This function requires that the inputs a and b be coprime.
PCG32_HOST_DEVICE static inline unsigned PCG32ModInverse(unsigned a,unsigned b);
// get nanoeconds from 1970.1.1 00:00:00 UTC
PCG32_HOST        static inline long long unsigned int PCG32TimeNanoeconds();
// 0 < x, no use
PCG32_HOST_DEVICE static inline double PCG32Ln(const double x);
PCG32_HOST_DEVICE static inline double PCG32PowerUnsigned(double x,unsigned n);
PCG32_HOST_DEVICE static inline double PCG32LnStirlingLeft(const double x);
// only accept x > 0
PCG32_HOST_DEVICE static inline double PCG32LnStirling(double x);


// Random Number Generator C Funtions (Prefix: PCG32)
PCG32_HOST_DEVICE static inline unsigned PCG32(PCG32Struct* status);
// Use PCG32SetSingleSeed for single seed, or PCG32SetMultipleSeeds for multiple seeds.
PCG32_HOST_DEVICE static inline void     PCG32SetSingleSeed(PCG32Struct* status,const long long unsigned int seed);
PCG32_HOST        static inline void     PCG32SetMultipleSeeds(PCG32Struct* statusArray,const long long unsigned int* baseSeeds,const unsigned count);
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform(PCG32Struct* status,unsigned min,unsigned max);
// max can not smaller than min and gap can not be power of 2
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_Strict(PCG32Struct* status,const unsigned min,const unsigned max);
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_StrictRangeUnchanged(PCG32Struct* status);
PCG32_HOST_DEVICE static inline void     PCG32UniformSetStrictRange(PCG32Struct* status,const unsigned min,const unsigned max);
// max can not smaller than min
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_MaxBiggerThanMin(PCG32Struct* status,const unsigned min,const unsigned max);
PCG32_HOST_DEVICE static inline double   PCG32UniformReal(PCG32Struct* status,const double min,const double max);
// undefined behavior if length of xy < 2
PCG32_HOST_DEVICE static inline void     PCG32RandomPointInCycle(PCG32Struct* status,const double radius,double* xy);
// undefined behavior if length of xy < 2
PCG32_HOST_DEVICE static inline void     PCG32RandomPointInCycle_NonCenter(PCG32Struct* status,const double radius,double* xy);
PCG32_HOST_DEVICE static inline double   PCG32StandardNormal(PCG32Struct* status);
// undefined behavior if length of random < 2
PCG32_HOST_DEVICE static inline void     PCG32StandardNormal2D(PCG32Struct* status,double random[2]);
// undefined behavior if length of random < 3
PCG32_HOST_DEVICE static inline void     PCG32StandardNormal3D(PCG32Struct* status,double random[3]);
// undefined behavior if length of random < n
PCG32_HOST_DEVICE static inline void     PCG32StandardNormalNDimension(PCG32Struct* status,const unsigned n,double* random);
// undefined behavior if length of xyz < 3 
PCG32_HOST_DEVICE static inline void     PCG32RandomPointInSphere3D(PCG32Struct* status,const double radius,double* xyz);
// undefined behavior if length of coordinate < dimension
PCG32_HOST_DEVICE static inline void     PCG32RandomPointInSphereNDimension(PCG32Struct* status,const double radius,const unsigned dimension,double* coordinate);
// Currently, only the algorithm for a >= 1 has been implemented
PCG32_HOST_DEVICE static inline bool     PCG32GammaInitialize(PCG32Struct* status,const double alpha,const double beta);
PCG32_HOST_DEVICE static inline double   PCG32Gamma(PCG32Struct* status);
PCG32_HOST_DEVICE static inline unsigned PCG32Binomial(PCG32Struct* status,double probability,const unsigned repeatUnsigned);
PCG32_HOST_DEVICE static inline double   PCG32Exponential(PCG32Struct* status,const double lambda);
PCG32_HOST_DEVICE static inline double   PCG32PowerLaw(PCG32Struct* status,const double min,const double alpha);
PCG32_HOST_DEVICE static inline bool     PCG32PoissonInitialize(PCG32Struct* status,const double mu);
PCG32_HOST_DEVICE static inline unsigned PCG32Poisson(PCG32Struct* status);
PCG32_HOST_DEVICE static inline double   PCG32StandardLogNormal(PCG32Struct* status);
PCG32_HOST_DEVICE static inline double   PCG32LogNormal(PCG32Struct* status,const double mu,const double sigma);
PCG32_HOST_DEVICE static inline double   PCG32Benford_SpecificLength(PCG32Struct* status,const unsigned minDigitalLength,const unsigned maxDigitalLength);
PCG32_HOST_DEVICE static inline double   PCG32Benford(PCG32Struct* status);
// n >= 1, undefined behivor if length of variables < n
PCG32_HOST_DEVICE static inline void     PCG32UniformSumReal(PCG32Struct* status,const unsigned n,const double sum,double* variables);
PCG32_HOST_DEVICE static inline unsigned PCG32Geometric(PCG32Struct* status,const double probability);
PCG32_HOST_DEVICE static inline unsigned PCG32Geometric_SmallProbability(PCG32Struct* status,const double probability);
// sum of probabilities must be 1
PCG32_HOST_DEVICE static inline unsigned PCG32MultinomialSampling(PCG32Struct* status,const double* probabilities,const unsigned length);
// sum of probabilities must be 1
PCG32_HOST_DEVICE static inline void     PCG32MultinomialSamplingCount(PCG32Struct* status,const double* probabilities,const unsigned length,const unsigned count,unsigned* result);


// AVX512 C Functions 
#if PCG32_AVX512

#include <immintrin.h>

// MUST BE BUILD WITH: g++/gcc -mavx512f -mavx512dq -mfma

typedef struct __x16__PCG32Struct{
    __m512i state0;
    __m512i state1;
    unsigned uniformStrictRange;
    unsigned uniformStrictM;
    unsigned uniformStrictShift;
    unsigned uniformStrictMin;
    unsigned uniformStrictGap;
}__attribute__((aligned(64))) __x16__PCG32Struct;

typedef __attribute__((aligned(64))) long long unsigned int __x16__StateArray[16];
typedef __attribute__((aligned(64))) long long unsigned int __x16__SeedArray[16];
typedef __attribute__((aligned(64))) unsigned               __x16__UnsignedArray[16];
typedef __attribute__((aligned(64))) double                 __x16__DoubleArray[16];

// Auxiliary C Functions
static inline __m512i U32VectorMultipleU32High32(const __m512i a,const unsigned b);
static inline __m512d __x8__PCG32Sin(const __m512d x);
static inline __m512d __x8__PCG32Cos(const __m512d x);
static inline __m512d __x8__PCG32Ln(const __m512d x);
static inline __m512d m512dRightShift(__m512d x,const unsigned shift64);

// Random Number Generator C Funtions (Prefix: __x16__PCG32)
static inline void __x16__PCG32(__x16__PCG32Struct* status,__x16__UnsignedArray random);
static inline void __x16__PCG32SetSingleSeed(__x16__PCG32Struct* status,__x16__SeedArray seed);
static inline void __x16__PCG32SetMultipleSeeds(__x16__PCG32Struct* status,__x16__SeedArray baseSeed[],const unsigned count);
static inline void __x16__PCG32UniformReal(__x16__PCG32Struct* status,__x16__DoubleArray random);
static inline void __x16__PCG32UniformReal_MinMax(__x16__PCG32Struct* status,const double min,const double max,__x16__DoubleArray random);
static inline void __x16__PCG32UniformSetStrictRange(__x16__PCG32Struct* status,const unsigned min,const unsigned max);
static inline void __x16__PCG32Uniform_StrictRangeUnchanged(__x16__PCG32Struct* status,__x16__UnsignedArray random);
static inline void __x16__PCG32Normal(__x16__PCG32Struct* status,__x16__DoubleArray random,const double mu,const double sigma);
static inline void __x16__PCG32StandardNormal(__x16__PCG32Struct* status,__x16__DoubleArray random);

#endif

#if defined(__cplusplus)||PCG32_CUDA

#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>

class PCG32PRNG{
public:
    PCG32_HOST_DEVICE PCG32PRNG();
    PCG32_HOST_DEVICE PCG32PRNG(long long unsigned int seed);
    PCG32_HOST_DEVICE ~PCG32PRNG(){}
    PCG32_HOST_DEVICE void SetSeed(long long unsigned int seed);

    PCG32_HOST_DEVICE unsigned Rand();
    PCG32_HOST_DEVICE unsigned Uniform(unsigned min,unsigned max);
    PCG32_HOST_DEVICE unsigned Uniform_Strict(const unsigned min,const unsigned max);
    PCG32_HOST_DEVICE void UniformSetStrictRange(const unsigned min,const unsigned max);
    PCG32_HOST_DEVICE unsigned Uniform_StrictRangeUnchanged();
    PCG32_HOST_DEVICE unsigned Uniform_MaxBiggerThanMin(const unsigned min,const unsigned max);
    PCG32_HOST_DEVICE double UniformReal(const double min,const double max);
    PCG32_HOST_DEVICE double Normal();
    PCG32_HOST_DEVICE double Normal(const double mu,const double sigma);
    PCG32_HOST_DEVICE bool GammaInitialize(const double alpha,const double beta);
    PCG32_HOST_DEVICE double Gamma();
    PCG32_HOST_DEVICE unsigned Binomial(double probability,const unsigned repeatUnsigned);
    template<typename Type>
    PCG32_HOST_DEVICE void UniformShuffle(Type* array,const long long unsigned int length);
    PCG32_HOST_DEVICE double Exponential(const double lambda);
    PCG32_HOST_DEVICE double PowerLaw(const double min,const double alpha);
    PCG32_HOST_DEVICE void RandomPointInCycle(const double radius,double* xy);
    PCG32_HOST_DEVICE double LogNormal();
    PCG32_HOST_DEVICE double LogNormal(const double mu,const double sigma);
    PCG32_HOST_DEVICE unsigned Geometric(const double probability);
    PCG32_HOST_DEVICE bool PoissonInitialize(const double mu);
    PCG32_HOST_DEVICE double Poisson();

    using result_type=uint32_t;
    PCG32_HOST_DEVICE static constexpr result_type min(){return 0;}
    PCG32_HOST_DEVICE static constexpr result_type max(){return PCG32MAX;}

    PCG32_HOST_DEVICE result_type operator()();
private:
    PCG32Struct status;
};

using PCG32PseudoRandomNumberGenerator=PCG32PRNG;

#endif

#ifdef __cplusplus
extern "C"{
#endif

PCG32_HOST_DEVICE static inline unsigned PCG32(PCG32Struct* status){
    long long unsigned int x=status->state;
    unsigned count=(unsigned)(x>>59);
    status->state=x*PCG32MULTIPLIER+PCG32INCREMENT;
    x=x^(x>>18);
    return rotr32((unsigned)(x>>27),count);
}

PCG32_DEPRECATED("Use PCG32SetSingleSeed for single seed, or PCG32SetMultipleSeeds for multiple seeds.")
PCG32_HOST_DEVICE static inline void PCG32SetSeed(PCG32Struct* status,const long long unsigned int seed){
    status->state=seed+PCG32INCREMENT;
    status->normalDistributionSavedValid=PCG32FALSE;
    PCG32(status);
}

PCG32_HOST_DEVICE static inline void PCG32SetSingleSeed(PCG32Struct* status,const long long unsigned int seed){
    status->state=seed+PCG32INCREMENT;
    status->normalDistributionSavedValid=PCG32FALSE;
    PCG32(status);
}

PCG32_HOST static inline void PCG32SetMultipleSeeds(PCG32Struct* statusArray,const long long unsigned int* baseSeeds,const unsigned count){
    unsigned prime=0;
    for(unsigned index=0;index<count;index=index+1){
        PCG32SetSingleSeed(statusArray+index,baseSeeds[index]);
        prime=PCG32NextPrime(prime);
        for(unsigned step=0;step<prime;step=step+1){
            PCG32(statusArray+index);
        }
    }
}

PCG32_HOST_DEVICE static inline unsigned PCG32Uniform(PCG32Struct* status,unsigned min,unsigned max){
    if(min>max){
        unsigned tempory=max;
        max=min;
        min=tempory;
    }
    long long unsigned int gap=max-min+1;
    if((gap&(gap-1))==0){
        return (PCG32(status)&(gap-1))+min;
    }
    unsigned range=(unsigned)(((PCG32MAX+1)/gap)*gap);
    unsigned random=PCG32(status);
    while(random>=range){
        random=PCG32(status);
    }
    return (random%gap)+min;
}

// max can not smaller than min and gap can not be power of 2
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_Strict(PCG32Struct* status,const unsigned min,const unsigned max){
    unsigned gap=max-min+1;
    unsigned range=(unsigned)(((PCG32MAX+1)/gap)*gap);
    unsigned random=PCG32(status);
    while(random>=range){
        random=PCG32(status);
    }
    return (random%gap)+min;
}

PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_StrictRangeUnchanged(PCG32Struct* status){
    unsigned random=PCG32(status);
    while(random>=status->uniformStrictRange){
        random=PCG32(status);
    }
    return (random%status->uniformStrictGap)+status->uniformStrictMin;
}

PCG32_HOST_DEVICE static inline void PCG32UniformSetStrictRange(PCG32Struct* status,const unsigned min,const unsigned max){
    unsigned gap=max-min+1;
    unsigned range=(unsigned)(((PCG32MAX+1)/gap)*gap);
    status->uniformStrictGap=gap;
    status->uniformStrictRange=range;
    status->uniformStrictMin=min;
}

// max can not smaller than min
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_MaxBiggerThanMin(PCG32Struct* status,const unsigned min,const unsigned max){
    long long unsigned int gap=max-min+1;
    if((gap&(gap-1))==0){
        return (PCG32(status)&(gap-1))+min;
    }
    unsigned range=(unsigned)(((PCG32MAX+1)/gap)*gap);
    unsigned random=PCG32(status);
    while(random>=range){
        random=PCG32(status);
    }
    return (random%gap)+min;
}

PCG32_HOST_DEVICE static inline double PCG32UniformReal(PCG32Struct* status,const double min,const double max){
    return min+((double)PCG32(status))*PCG32REAL_SCALE*(max-min);
}

// undefined behavior if length of xy < 2
PCG32_HOST_DEVICE static inline void PCG32RandomPointInCycle(PCG32Struct* status,const double radius,double* xy){
    const double Pi2=6.2831853071795864769252867665590057683943387987502L;
    const double alpha=PCG32UniformReal(status,0,Pi2);
    const double r=radius*sqrt(PCG32UniformReal(status,0,1));
    xy[0]=r*cos(alpha);
    xy[1]=r*sin(alpha);
}

// undefined behavior if length of xy < 2
PCG32_HOST_DEVICE static inline void PCG32RandomPointInCycle_NonCenter(PCG32Struct* status,const double radius,double* xy){
    const double Pi2=6.2831853071795864769252867665590057683943387987502L;
    const double alpha=PCG32UniformReal(status,0,Pi2);
    const double r=radius*sqrt(1-PCG32UniformReal(status,0,1));
    xy[0]=r*cos(alpha);
    xy[1]=r*sin(alpha);
}

PCG32_HOST_DEVICE static inline double PCG32StandardNormal(PCG32Struct* status){
    if(status->normalDistributionSavedValid){
        status->normalDistributionSavedValid=PCG32FALSE;
        return status->normalDistributionSaved;
    }
    #if 1
    double u1,u2,S;
    do{
        u1=PCG32UniformReal(status,-1,1);
        u2=PCG32UniformReal(status,-1,1);
        S=u1*u1+u2*u2;
    }while(S>1.0||S==0.0);
    #else
    double u[2];
    PCG32RandomPointInCycle_NonCenter(status,1,u);
    const double u1=u[0];
    const double u2=u[1];
    const double S=u1*u1+u2*u2;
    #endif
    const double toMultiple=sqrt(-2.0*log(S)/S);
    status->normalDistributionSavedValid=PCG32TRUE;
    status->normalDistributionSaved=toMultiple*u2;
    return toMultiple*u1;
}

// undefined behavior if length of random < 2
PCG32_HOST_DEVICE static inline void PCG32StandardNormal2D(PCG32Struct* status,double random[2]){
    double u1,u2,S;
    do{
        u1=PCG32UniformReal(status,-1,1);
        u2=PCG32UniformReal(status,-1,1);
        S=u1*u1+u2*u2;
    }while(S>1.0||S==0.0);
    const double toMultiple=sqrt(-2.0*log(S)/S);
    random[0]=toMultiple*u1;
    random[1]=toMultiple*u2;
}

// undefined behavior if length of random < 3
PCG32_HOST_DEVICE static inline void PCG32StandardNormal3D(PCG32Struct* status,double random[3]){
    double u1,u2,S;
    do{
        u1=PCG32UniformReal(status,-1,1);
        u2=PCG32UniformReal(status,-1,1);
        S=u1*u1+u2*u2;
    }while(S>1.0||S==0.0);
    const double toMultiple=sqrt(-2.0*log(S)/S);
    random[0]=toMultiple*u1;
    random[1]=toMultiple*u2;
    random[2]=PCG32StandardNormal(status);
}

// undefined behavior if length of random < n
PCG32_HOST_DEVICE static inline void PCG32StandardNormalNDimension(PCG32Struct* status,const unsigned n,double* random){
    for(unsigned index=0;index<n;index=index+2){
        double u1,u2,S;
        do{
            u1=PCG32UniformReal(status,-1,1);
            u2=PCG32UniformReal(status,-1,1);
            S=u1*u1+u2*u2;
        }while(S>1.0||S==0.0);
        const double toMultiple=sqrt(-2.0*log(S)/S);
        random[index+0]=toMultiple*u1;
        random[index+1]=toMultiple*u2;
    }
    if(n&1){
        random[n-1]=PCG32StandardNormal(status);
    }
}

// undefined behavior if length of xyz < 3 
PCG32_HOST_DEVICE static inline void PCG32RandomPointInSphere3D(PCG32Struct* status,const double radius,double* xyz){
    xyz[0]=PCG32StandardNormal(status);
    xyz[1]=PCG32StandardNormal(status);
    xyz[2]=PCG32StandardNormal(status);
    const double length=1.0/sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
    const double r=radius*pow(PCG32UniformReal(status,0,1),1.0/3.0);
    xyz[0]=xyz[0]*r*length;
    xyz[1]=xyz[1]*r*length;
    xyz[2]=xyz[2]*r*length;
}

// undefined behavior if length of coordinate < dimension
PCG32_HOST_DEVICE static inline void PCG32RandomPointInSphereNDimension(PCG32Struct* status,const double radius,const unsigned dimension,double* coordinate){
    double sum=0.0;
    for(unsigned index=0;index<dimension;index=index+1){
        coordinate[index]=PCG32StandardNormal(status);
        sum=sum+coordinate[index]*coordinate[index];
    }
    const double length=1.0/sqrt(sum);
    const double r=radius*pow(PCG32UniformReal(status,0,1),1.0/dimension);
    for(unsigned index=0;index<dimension;index=index+1){
        coordinate[index]=coordinate[index]*r*length;
    }
}

// Currently, only the algorithm for a >= 1 has been implemented
PCG32_HOST_DEVICE static inline bool PCG32GammaInitialize(PCG32Struct* status,const double alpha,const double beta){
    if(alpha<1.0){
        return false;
    }
    status->gammaD=alpha-1.0/3.0;
    status->gammaC=1.0/sqrt(9*status->gammaD);
    status->gammaNegativeSqrt9AlphaSub3=-1.0*sqrt(9.0*alpha-3.0);
    if(beta>0.0){
        status->gammaBeta=beta;
        return true;
    }
    return false;
}

PCG32_HOST_DEVICE static inline double PCG32Gamma(PCG32Struct* status){
    double normal,normal2,normal4;
    const double NegativeSqrt9AlphaSub3=status->gammaNegativeSqrt9AlphaSub3;
    while(true){
        do{
            normal=PCG32StandardNormal(status);
        }while(normal<NegativeSqrt9AlphaSub3);
        double v=(1+status->gammaC*normal);
        v=v*v*v;
        double uniform=PCG32UniformReal(status,0,1);
        normal2=normal*normal;
        normal4=1-normal2*normal2*0.0331;
        if(uniform<normal4){
            return status->gammaD*v*status->gammaBeta;
        }
        if(log(uniform)<(normal2*0.5+status->gammaD*(1-v+log(v)))){
            return status->gammaD*v*status->gammaBeta;
        }
    }
}

PCG32_HOST_DEVICE static inline unsigned PCG32Binomial(PCG32Struct* status,double probability,const unsigned repeatUnsigned){
    const unsigned PCG32BINOMIAL_SMALLMEAN=14;
    const unsigned PCG32BINOMIAL_MAXITERATION=110;
    const unsigned PCG32BINOMIAL_FARFROMMENA=20;
    int repeats[2]={0,0};
    unsigned result=0;
    if(repeatUnsigned>PCG32_INT_MAX){
        if(repeatUnsigned==PCG32_UNSIGNED_MAX){
            repeats[0]=PCG32_INT_MAX;
            repeats[1]=PCG32_INT_MAX;
            for(unsigned index=0;index<(PCG32_UNSIGNED_MAX-PCG32_INT_MAX-PCG32_INT_MAX);index=index+1){
                if(PCG32UniformReal(status,0,1)<probability){
                    result=result+1;
                }
            }
        }else{
            repeats[0]=PCG32_INT_MAX;
            repeats[1]=repeatUnsigned-PCG32_INT_MAX;
        }
    }else{
        repeats[0]=repeatUnsigned;
    }
    for(unsigned repeatsIndex=0;repeatsIndex<2;repeatsIndex=repeatsIndex+1){
        unsigned repeat=repeats[repeatsIndex];
        if(!repeat){
            break;
        }
        int success=0;
        bool flipped=false;
        if(probability>0.5){
            probability=1.0-probability;
            flipped=true;
        }
        const double q=1.0-probability;
        const double rate=probability/q;
        const double mean=repeat*probability;
        if(mean<PCG32BINOMIAL_SMALLMEAN){
            double f0=PCG32PowerUnsigned(q,repeat);
            while(true){
                double f=f0;
                double uniform=PCG32UniformReal(status,0,1);
                for(success=0;success<=PCG32BINOMIAL_MAXITERATION;success=success+1){
                    if(uniform<f){
                        goto Finish;
                    }
                    uniform=uniform-f;
                    f=f*rate*((double)(repeat-success)/(double)(success+1));
                }
            }
        }else{
            int k;
            const double ffm=mean+probability;
            const int m=(int)ffm;
            const double fm=m;
            const double xm=fm+0.5;
            const double meanq=mean*q;
            const double c=0.134+20.5/(15.3+fm);
            double p[4];
            p[0]=floor(2.195*sqrt(meanq)-4.6*q)+0.5;
            p[1]=p[0]*(1.0+c+c);
            const double xl=xm-p[0];
            const double xr=xm+p[0];
            const double al=(ffm-xl)/(ffm-xl*probability);
            const double ar=(xr-ffm)/(xr*q);
            const double lambdaL=al*(1.0+0.5*al);
            const double lambdaR=ar*(1.0+0.5*ar);
            p[2]=p[1]+c/lambdaL;
            p[3]=p[2]+c/lambdaR;

            double var,accept,u,v;
            TryAgain:
            u=PCG32UniformReal(status,0,p[3]);
            v=PCG32UniformReal(status,0,1);
            if(u<=p[0]){
                success=(int)(xm-p[0]*v+u);
                goto Finish;
            }else if(u<=p[1]){
                const double x=xl+(u-p[0])/c;
                v=v*c+1.0-fabs(x-xm)/p[0];
                if(v>1.0||v<=0.0){
                    goto TryAgain;
                }
                success=(int)x;
            }else if(u<=p[2]){
                int successInt=(int)(xl+log(v)/lambdaL);
                if(successInt<0){
                    goto TryAgain;
                }
                success=successInt;
                v=v*((u-p[1])*lambdaL);
            }else{
                success=(int)(xr-log(v)/lambdaR);
                if(success>repeat){
                    goto TryAgain;
                }
                v=v*((u-p[2])*lambdaR);
            }
            if(success>m){
                k=success-m;
            }else{
                k=m-success;
            }
            if(k<=PCG32BINOMIAL_FARFROMMENA){
                const double g=(repeat+1)*rate;
                double f=1.0;
                var=v;
                if(m<success){
                    for(int index=m+1;index<=success;index=index+1){
                        f=f*(g/index-rate);
                    }
                }else if(m>success){
                    for(int index=success+1;index<=m;index=index+1){
                        f=f/(g/index-rate);
                    }
                }
                accept=f;
            }else{
                var=log(v);
                if(k<meanq/2-1){
                    const double amaxp=k/meanq*((k*(k/3.0+0.625)+(1.0/6.0))/meanq+0.5);
                    const double ynorm=-((k*k)/(2.0*meanq));
                    if(var<ynorm-amaxp){
                        goto Finish;
                    }
                    if(var>ynorm+amaxp){
                        goto TryAgain;
                    }
                }
                const double x1=success+1.0;
                const double w1=repeat-success+1.0;
                const double f1=fm+1.0;
                const double z1=repeat+1.0-fm;
                accept=xm*log(f1/x1)+(repeat-m+0.5)*log(z1/w1)+(success-m)*log(w1*probability/(x1*q))+PCG32LnStirlingLeft(f1)+PCG32LnStirlingLeft(z1)-PCG32LnStirlingLeft(x1)-PCG32LnStirlingLeft(w1);
            }
            if(var<=accept){
                goto Finish;
            }else{
                goto TryAgain;
            }
        }
        Finish:
        if(flipped){
            result=result+(repeat-success);
        }else{
            result=result+success;
        }
    }
    return result;
}

PCG32_HOST_DEVICE static inline double PCG32Exponential(PCG32Struct* status,const double lambda){
    return -log(1-PCG32UniformReal(status,0,1))/lambda;
}

PCG32_HOST_DEVICE static inline double PCG32PowerLaw(PCG32Struct* status,const double min,const double alpha){
    return min*pow(1-PCG32UniformReal(status,0,1),1.0/(1.0-alpha));
}

PCG32_HOST_DEVICE static inline bool PCG32PoissonInitialize(PCG32Struct* status,const double mu){
    if(mu<=0.0){
        return false;
    }
    if(mu>=12){
        const double floorMu=floor(mu);
        const double PiDive4  =0.7853981633974483096156608458198757210492923498437764L;
        const double PiDiveBy4=1.2732395447351626861510701069801148962756771659236516L;
        const double dx=sqrt(2*floorMu*log(32*floorMu*PiDiveBy4));
        #if defined(__cplusplus)&&(!PCG32_CUDA)
        constexpr double epsilon=(1.0-std::numeric_limits<double>::epsilon())/2;
        const double roundMu=std::round(std::max(6.0,std::min(floorMu,dx)));
        #else
        double roundMuTemporary;
        double min;
        if(floorMu<dx){
            min=floorMu;
        }else{
            min=dx;
        }
        if(6.0<min){
            roundMuTemporary=round(min);
        }else{
            roundMuTemporary=6.0;
        }
        const double roundMu=roundMuTemporary;
        #endif
        const double sqrtPi2 =1.2533141373155002512078826424055226265034933703049692L;
        const double logMu=log(mu);
        const double lgmma=PCG32LnStirling(floorMu);
        const double sqrtMu=sqrt(floorMu);
        double c[6];
        c[1]=sqrtMu*sqrtPi2;
        const double cx=2*floorMu+roundMu;
        const double sqrtCx=sqrt(cx/2);
        const double cxDiviedBy1=1.0/cx;
        const double c2b=sqrt(PiDive4*cx)*exp(cxDiviedBy1);
        const double cb=2*cx*exp(-roundMu*cxDiviedBy1*(1.0+roundMu/2.0))/roundMu;
        c[2]=c2b+c[1];
        c[3]=c[2]+1;
        c[4]=c[3]+1;
        const double _1dive78=0.0128205128205128205128205128205128205128205128205128L;
        const double e1dive78=1.0129030479320018583185514777512982888868114426296881L;
        c[5]=c[4]+e1dive78;
        c[0]=cb+c[5];
        const double cx2=4*floorMu+2*roundMu;

        status->Poisson.mu=mu;
        status->Poisson.floorMu=floorMu;
        status->Poisson.dx=dx;
        status->Poisson.roundMu=roundMu;
        status->Poisson.logMu=logMu;
        status->Poisson.lgmma=lgmma;
        status->Poisson.sqrtMu=sqrtMu;
        for(unsigned index=0;index<6;index=index+1){
            status->Poisson.c[index]=c[index];
        }
        status->Poisson.cx=cx;
        status->Poisson.sqrtCx=sqrtCx;
        status->Poisson.cxDiviedBy1=cxDiviedBy1;
        status->Poisson.c2b=c2b;
        status->Poisson.cb=cb;
        status->Poisson.cx2=cx2;
    }else{
        const double expmu=exp(-mu);

        status->Poisson.mu=mu;
        status->Poisson.expmu=expmu;
    }
    return true;
}

PCG32_HOST_DEVICE static inline unsigned PCG32Poisson(PCG32Struct* status){
    const double mu=status->Poisson.mu;
    if(mu>=12){
        double k;
        const double PiDive4  =0.7853981633974483096156608458198757210492923498437764L;
        const double PiDiveBy4=1.2732395447351626861510701069801148962756771659236516L;
        #if defined(__cplusplus)&&(!PCG32_CUDA)
        constexpr double epsilon=(1.0-std::numeric_limits<double>::epsilon())/2;
        constexpr double threshod=(double)PCG32MAX+epsilon;
        #else
        const double epsilon=(1.0-DBL_EPSILON)/2;
        const double threshod=(double)PCG32MAX+epsilon;
        #endif
        const double sqrtPi2 =1.2533141373155002512078826424055226265034933703049692L;
        const double _1dive78=0.0128205128205128205128205128205128205128205128205128L;
        const double e1dive78=1.0129030479320018583185514777512982888868114426296881L;

        const double floorMu=status->Poisson.floorMu;
        const double dx=status->Poisson.dx;
        const double roundMu=status->Poisson.roundMu;
        const double logMu=status->Poisson.logMu;
        const double lgmma=status->Poisson.lgmma;
        const double sqrtMu=status->Poisson.sqrtMu;
        const double c[]={status->Poisson.c[0],status->Poisson.c[1],status->Poisson.c[2],status->Poisson.c[3],status->Poisson.c[4],status->Poisson.c[5]};
        const double cx=status->Poisson.cx;
        const double sqrtCx=status->Poisson.sqrtCx;
        const double cxDiviedBy1=status->Poisson.cxDiviedBy1;
        const double c2b=status->Poisson.c2b;
        const double cb=status->Poisson.cb;
        const double cx2=status->Poisson.cx2;

        while(true){
            const double uniform=c[0]*PCG32UniformReal(status,0,1);
            const double ln=-log(1.0-PCG32UniformReal(status,0,1));
            double w=0.0;
            if(uniform<=c[1]){
                const double normal=PCG32StandardNormal(status);
                k=floor(-fabs(normal)*sqrtMu-1);
                w=normal/2*(-normal);
                if(k<-floorMu){
                    continue;
                }
            }else if(uniform<=c[2]){
                const double normal=PCG32StandardNormal(status);
                const double unceil=fabs(normal)*sqrtCx+1;
                k=ceil(unceil);
                w=unceil*(2-unceil)*cxDiviedBy1;
                if(k>roundMu){
                    continue;
                }
            }else if(uniform<=c[3]){
                k=-1;
            }else if(uniform<=c[4]){
                k=0;
            }else if(uniform<=c[5]){
                k=1;
                w=_1dive78;
            }else{
                const double lnelse=-log(1.0-PCG32UniformReal(status,0,1));
                const double unceil=roundMu+lnelse*cx2/roundMu;
                k=ceil(unceil);
                w=-roundMu*cxDiviedBy1*(1+unceil/2);
            }
            if(w-ln-k*logMu<=lgmma-PCG32LnStirling(k+floorMu)){
                if(k+floorMu<threshod){
                    break;
                }
            }
        }
        return (unsigned)(k+floorMu+epsilon);
    }else{
        unsigned k=0;
        double product=1.0;
        const double expmu=status->Poisson.expmu;
        while(true){
            product=product*PCG32UniformReal(status,0,1);
            if(product<=expmu){
                return k;
            }
            k=k+1;
        }
    }
}

PCG32_HOST_DEVICE static inline double PCG32StandardLogNormal(PCG32Struct* status){
    return exp(PCG32StandardNormal(status));
}

PCG32_HOST_DEVICE static inline double PCG32LogNormal(PCG32Struct* status,const double mu,const double sigma){
    return exp(mu+PCG32StandardNormal(status)*sigma);
}

PCG32_HOST_DEVICE static inline double PCG32Benford_SpecificLength(PCG32Struct* status,const unsigned minDigitalLength,const unsigned maxDigitalLength){
    return pow(10,PCG32UniformReal(status,minDigitalLength-1,maxDigitalLength));
}

PCG32_HOST_DEVICE static inline double PCG32Benford(PCG32Struct* status){
    return PCG32Benford_SpecificLength(status,1,12);
}

// n >= 1
// n uniform real variables that with sum to sum
// undefined behivor if length of variables < n
PCG32_HOST_DEVICE static inline void PCG32UniformSumReal(PCG32Struct* status,const unsigned n,const double sum,double* variables){
    #if 0
    double R=sum;
    for(unsigned index=0;index<n-1;index=index+1){
        const double nextR=R*pow(PCG32UniformReal(status,0,1),1.0/(n-index-1));
        variables[index]=R-nextR;
        R=nextR;
    }
    variables[n-1]=R;
    #else
    double S=0;
    for(unsigned index=0;index<n;index=index+1){
        variables[index]=PCG32Exponential(status,1);
        S=S+variables[index];
    }
    const double scale=sum/S;
    for(unsigned index=0;index<n;index=index+1){
        variables[index]=variables[index]*scale;
    }
    #endif
}

PCG32_HOST_DEVICE static inline unsigned PCG32Geometric(PCG32Struct* status,const double probability){
    const double uniform=PCG32UniformReal(status,0,1);
    return (unsigned)(log(uniform)/log(1.0-probability))+1;
}

PCG32_HOST_DEVICE static inline unsigned PCG32Geometric_SmallProbability(PCG32Struct* status,const double probability){
    const double uniform=PCG32UniformReal(status,0,1);
    return (unsigned)(log(uniform)/log1p(-probability))+1;
}

// sum of probabilities must be 1
PCG32_HOST_DEVICE static inline unsigned PCG32MultinomialSampling(PCG32Struct* status,const double* probabilities,const unsigned length){
    const double uniform=PCG32UniformReal(status,0,1);
    if(probabilities[0]>=uniform){
        return 0;
    }
    double sum=probabilities[0];
    for(unsigned index=1;index<length;index=index+1){
        sum=sum+probabilities[index];
        if(sum>=uniform){
            return index;
        }
    }
}

// sum of probabilities must be 1
PCG32_HOST_DEVICE static inline void PCG32MultinomialSamplingCount(PCG32Struct* status,const double* probabilities,const unsigned length,const unsigned count,unsigned* result){
    double probabilityLeft=1;
    unsigned sumCount=0;
    for(unsigned index=0;index<length-1;index=index+1){
        result[index]=PCG32Binomial(status,probabilities[index]/probabilityLeft,count-sumCount);
        sumCount=sumCount+result[index];
        probabilityLeft=probabilityLeft-probabilities[index];
    }
    result[length-1]=count-sumCount;
}

#if PCG32_AVX512

#include <immintrin.h>
#include <stdatomic.h> // for m512dLeftShift compile by GCC 

// MUST BE BUILD WITH: g++/gcc -mavx512f -mavx512dq -mfma

#define PCG32_AVX512_PCG32_INTRINSIC_LOAD_CONSTANT \
    const __m512i add     =_mm512_set1_epi64(PCG32INCREMENT);                                 \
    const __m512i multiple=_mm512_set1_epi64(PCG32MULTIPLIER);                                \
    const __m512i toAnd   =_mm512_set1_epi64(0x1F0000001F);                                   

#define PCG32_AVX512_PCG32_INTRINSIC_CALCULATE \
    __m512i x;                                                                                \
    {                                                                                         \
        __m512i x0=status->state0;                                                            \
        __m512i x1=status->state1;                                                            \
        const __m512i count0=_mm512_srli_epi64(x0,59);                                        \
        const __m512i count1=_mm512_srli_epi64(x1,59);                                        \
        status->state0=_mm512_add_epi64(_mm512_mullo_epi64(x0,multiple),add);                 \
        status->state1=_mm512_add_epi64(_mm512_mullo_epi64(x1,multiple),add);                 \
        x0=x0^_mm512_srli_epi64(x0,18);                                                       \
        x1=x1^_mm512_srli_epi64(x1,18);                                                       \
        x0=_mm512_srli_epi64(x0,27);                                                          \
        x1=_mm512_srli_epi64(x1,27);                                                          \
        const __m256i count0Unsigned=_mm512_cvtepi64_epi32(count0);                           \
        const __m256i count1Unsigned=_mm512_cvtepi64_epi32(count1);                           \
        __m512i r;                                                                            \
        r=_mm512_inserti64x4(r,count0Unsigned,0);                                             \
        r=_mm512_inserti64x4(r,count1Unsigned,1);                                             \
        const __m256i x0Unsigned=_mm512_cvtepi64_epi32(x0);                                   \
        const __m256i x1Unsigned=_mm512_cvtepi64_epi32(x1);                                   \
        x=_mm512_inserti64x4(x,x0Unsigned,0);                                                 \
        x=_mm512_inserti64x4(x,x1Unsigned,1);                                                 \
        const __m512i zero=_mm512_setzero_epi32();                                            \
        const __m512i _r=_mm512_sub_epi32(zero,r);                                            \
        const __m512i right=_mm512_srlv_epi32(x,r);                                           \
        const __m512i left =_mm512_sllv_epi32(x,_mm512_and_epi32(_r,toAnd));                  \
        x=_mm512_or_epi32(right,left);                                                        \
    }

#define PCG32_AVX512_PCG32_INTRINSIC \
    PCG32_AVX512_PCG32_INTRINSIC_LOAD_CONSTANT \
    PCG32_AVX512_PCG32_INTRINSIC_CALCULATE

#define PCG32_AVX512_UNIFORM_REAL_INTRINSIC \
    __m512d random0,random1;                                                                  \
    {                                                                                         \
        const __m256i x0=_mm512_extracti64x4_epi64(x,0);                                      \
        const __m256i x1=_mm512_extracti64x4_epi64(x,1);                                      \
        const __m512d x0Double=_mm512_cvtepu32_pd(x0);                                        \
        const __m512d x1Double=_mm512_cvtepu32_pd(x1);                                        \
        const __m512i multipleDouble=_mm512_set1_epi64(0x3DF0000000000000LLU);                \
        random0=_mm512_mul_pd(x0Double,_mm512_castsi512_pd(multipleDouble));                  \
        random1=_mm512_mul_pd(x1Double,_mm512_castsi512_pd(multipleDouble));                  \
    }

#define PCG32_AVX512_UNIFORM_REAL_MIN_MAX_INTRINSIC(min,max) \
    {                                                                                         \
        union DoubleLLU{                                                                      \
            double d;                                                                         \
            long long unsigned int u;                                                         \
        }MinUnion,MaxUnion;                                                                   \
        MinUnion.d=min;                                                                       \
        MaxUnion.d=max;                                                                       \
        const __m512d minv=_mm512_castsi512_pd(_mm512_set1_epi64(MinUnion.u));                \
        const __m512d maxv=_mm512_castsi512_pd(_mm512_set1_epi64(MaxUnion.u));                \
        const __m512d gapv=_mm512_sub_pd(maxv,minv);                                          \
        random0=_mm512_mul_pd(random0,gapv);                                                  \
        random1=_mm512_mul_pd(random1,gapv);                                                  \
        random0=_mm512_add_pd(random0,minv);                                                  \
        random1=_mm512_add_pd(random1,minv);                                                  \
    }

static inline void __x16__PCG32(__x16__PCG32Struct* status,__x16__UnsignedArray random){
    PCG32_AVX512_PCG32_INTRINSIC
    _mm512_store_epi64(random,x);
}

static inline void __x16__PCG32SetSingleSeed(__x16__PCG32Struct* status,__x16__SeedArray seed){
    const __m512i seed0=_mm512_set_epi64(seed[ 7],seed[ 6],seed[ 5],seed[ 4],seed[ 3],seed[ 2],seed[ 1],seed[ 0]);
    const __m512i seed1=_mm512_set_epi64(seed[15],seed[14],seed[13],seed[12],seed[11],seed[10],seed[ 9],seed[ 8]);
    const __m512i add=_mm512_set_epi64(PCG32INCREMENT,PCG32INCREMENT,PCG32INCREMENT,PCG32INCREMENT,PCG32INCREMENT,PCG32INCREMENT,PCG32INCREMENT,PCG32INCREMENT);
    status->state0=_mm512_add_epi64(seed0,add);
    status->state1=_mm512_add_epi64(seed1,add);
    __x16__UnsignedArray random;
    __x16__PCG32(status,random);
}

static inline void __x16__PCG32SetMultipleSeeds(__x16__PCG32Struct* status,__x16__SeedArray baseSeed[],const unsigned count){
    PCG32Struct scalarStatus[16];
    unsigned prime=0;
    for(unsigned indexStatus=0;indexStatus<count;indexStatus=indexStatus+1){
        for(unsigned index=0;index<sizeof(__x16__SeedArray)/sizeof(baseSeed[0][0]);index=index+1){
            PCG32SetSingleSeed(scalarStatus+index,baseSeed[indexStatus][index]);
            prime=PCG32NextPrime(prime);
            for(unsigned step=0;step<prime;step=step+1){
                PCG32(scalarStatus+index);
            }
            baseSeed[indexStatus][index]=scalarStatus[index].state;
        }
        __x16__PCG32SetSingleSeed(status+indexStatus,baseSeed[indexStatus]);
    }
}

// 0 <= random < 1
static inline void __x16__PCG32UniformReal(__x16__PCG32Struct* status,__x16__DoubleArray random){
    PCG32_AVX512_PCG32_INTRINSIC
    PCG32_AVX512_UNIFORM_REAL_INTRINSIC
    _mm512_store_pd(random+0,random0);
    _mm512_store_pd(random+8,random1);
}

static inline void __x16__PCG32UniformReal_MinMax(__x16__PCG32Struct* status,const double min,const double max,__x16__DoubleArray random){
    PCG32_AVX512_PCG32_INTRINSIC                                                              
    PCG32_AVX512_UNIFORM_REAL_INTRINSIC                                                
    PCG32_AVX512_UNIFORM_REAL_MIN_MAX_INTRINSIC(min,max)
    _mm512_store_pd(random+0,random0);
    _mm512_store_pd(random+8,random1);
}

static inline __m512i U32VectorMultipleU32High32(const __m512i a,const unsigned b){
    const __m512i bv=_mm512_set1_epi32(b);
    const __m512i high0=_mm512_srli_epi64(_mm512_mul_epu32(a,bv),32);
    const __m512i aShift32=_mm512_srli_epi64(a,32);
    const __m512i mask=_mm512_set1_epi64(0xFFFFFFFF00000000LLU);
    const __m512i high1=_mm512_and_si512(_mm512_mul_epu32(aShift32,bv),mask);
    return _mm512_or_si512(high0,high1);
}

#define PCG32_AVX512_UNIFORM_STRICT_RANGE_UNCHANGED_MOD_METHOD 1

static inline void __x16__PCG32Uniform_StrictRangeUnchanged(__x16__PCG32Struct* status,__x16__UnsignedArray random){
    #if 1
    PCG32_AVX512_PCG32_INTRINSIC_LOAD_CONSTANT
    __m512i random16;
    const __m512i range=_mm512_set1_epi32(status->uniformStrictRange);
    bool first=true;
    while(true){
        PCG32_AVX512_PCG32_INTRINSIC_CALCULATE
        const __mmask16 compare=_mm512_cmplt_epu32_mask(x,range);
        const unsigned count=__builtin_popcount(compare);
        if(count>=8){
            if(count==16){
                random16=x;
                break;
            }
            if(first){
                random16=_mm512_mask_compress_epi32(_mm512_setzero_si512(),compare,x);
            }else{
                random16=_mm512_inserti64x4(random16,_mm512_castsi512_si256(_mm512_mask_compress_epi32(_mm512_setzero_si512(),compare,x)),1);
                break;
            }
            first=false;
        }
    }
    #else
    PCG32_AVX512_PCG32_INTRINSIC
    __m512i random16=x;
    #endif
    #if PCG32_AVX512_UNIFORM_STRICT_RANGE_UNCHANGED_MOD_METHOD==0

    unsigned shiftCount=status->shiftCount;
    __m512i gapShift=status->uniformStrictGapShift;
    while(shiftCount){
        const __mmask16 mask=_mm512_cmpge_epu32_mask(random16,gapShift);
        random16=_mm512_mask_sub_epi32(random16,mask,random16,gapShift);
        gapShift=_mm512_srli_epi32(gapShift,1);
        shiftCount=shiftCount-1;
    }
    const __m512i min=status->uniformStrictMin;
    random16=_mm512_add_epi32(random16,min);
    _mm512_store_epi64(random,random16);

    #elif PCG32_AVX512_UNIFORM_STRICT_RANGE_UNCHANGED_MOD_METHOD==1

    __m512i quotient=U32VectorMultipleU32High32(random16,status->uniformStrictM);
    quotient=_mm512_add_epi32(_mm512_srli_epi32(_mm512_sub_epi32(random16,quotient),1),quotient);
    quotient=_mm512_srli_epi32(quotient,status->uniformStrictShift);
    const __m512i gap=_mm512_set1_epi32(status->uniformStrictGap);
    const __m512i integerMultiple=_mm512_mullo_epi32(quotient,gap);
    random16=_mm512_sub_epi32(random16,integerMultiple);
    const __m512i min=_mm512_set1_epi32(status->uniformStrictMin);
    random16=_mm512_add_epi32(random16,min);
    _mm512_store_epi64(random,random16);

    #elif PCG32_AVX512_UNIFORM_STRICT_RANGE_UNCHANGED_MOD_METHOD==2

    _mm512_store_epi64(random,random16);
    const unsigned min=status->uniformStrictMin;
    const unsigned gap=status->uniformStrictGap;
    random[ 0]=random[ 0]%gap+min;
    random[ 1]=random[ 1]%gap+min;
    random[ 2]=random[ 2]%gap+min;
    random[ 3]=random[ 3]%gap+min;
    random[ 4]=random[ 4]%gap+min;
    random[ 5]=random[ 5]%gap+min;
    random[ 6]=random[ 6]%gap+min;
    random[ 7]=random[ 7]%gap+min;
    random[ 8]=random[ 8]%gap+min;
    random[ 9]=random[ 9]%gap+min;
    random[10]=random[10]%gap+min;
    random[11]=random[11]%gap+min;
    random[12]=random[12]%gap+min;
    random[13]=random[13]%gap+min;
    random[14]=random[14]%gap+min;
    random[15]=random[15]%gap+min;

    #elif PCG32_AVX512_UNIFORM_STRICT_RANGE_UNCHANGED_MOD_METHOD==3

    const __m512i modA=_mm512_set1_epi32((1<<status->uniformStrictMoreShift)-1);
    const __m512i r1=_mm512_and_si512(random16,modA);
    const __m512i min=_mm512_add_epi32(r1,_mm512_set1_epi32(status->uniformStrictMin));
    __m512i quotient=U32VectorMultipleU32High32(random16,status->uniformStrictM);
    quotient=_mm512_srli_epi32(quotient,status->uniformStrictShift);
    const __m512i gap=_mm512_set1_epi32(status->uniformStrictGap);
    const __m512i integerMultiple=_mm512_mullo_epi32(quotient,gap);
    random16=_mm512_sub_epi32(random16,integerMultiple);
    random16=_mm512_sub_epi32(random16,r1);
    const __m512i inverse=_mm512_set1_epi32(status->uniformStrictInverse);
    random16=_mm512_mullo_epi32(random16,inverse);
    random16=_mm512_slli_epi32(random16,status->uniformStrictMoreShift);
    random16=_mm512_add_epi32(random16,min);
    _mm512_store_epi64(random,r1);

    #else
        #error "Unspecified mod method. "
    #endif
}

static inline void __x16__PCG32UniformSetStrictRange(__x16__PCG32Struct* status,const unsigned min,const unsigned max){
    const unsigned gap=max-min+1;
    const unsigned range=(unsigned)(((PCG32MAX+1)/gap)*gap);
    unsigned shift;
    unsigned floorLog2Gap=31-__builtin_clz(gap);
    long long unsigned int M=((1LLU<<floorLog2Gap)<<32)/(long long unsigned int)gap;
    unsigned left=(unsigned)(((1LLU<<floorLog2Gap)<<32)-((long long unsigned int)gap*M));
    M=M<<1;
    const unsigned leftShift1=left<<1;
    if(leftShift1>=gap||leftShift1<left){
        M=M+1;
    }
    shift=floorLog2Gap;
    M=M+1;
    status->uniformStrictRange=range;
    status->uniformStrictShift=(unsigned)shift;
    status->uniformStrictMin=min;
    status->uniformStrictM=(unsigned)M;
    status->uniformStrictGap=gap;
}

#define PCG32_AVX512_SET_DOUBLE_LLU(x) _mm512_castsi512_pd(_mm512_set1_epi64(x))

/**
 * -pi/2 <= x <= pi/2
 * x1*(0.99999999999999999974277490079943975 + 
 * x2*(-0.166666666666666650522767323353840604 + 
 * x2*(0.00833333333333316503140948668861163462 + 
 * x2*(-0.00019841269841201840459252750531485886 + 
 * x2*(2.75573192101527564362114785169078252e-6 + 
 * x2*(-2.50521067982746148969440582709985054e-8 + 
 * x2*(1.60589364903732230834314189302038183e-10 + 
 * x2*(-7.64291780693694318128770390349958602e-13 + 
 * 2.72047909631134875287705126898888084e-15*x2))))))))
 */
static inline __m512d __x8__PCG32Sin(const __m512d x){
    const __m512d x1=x;
    const __m512d x2=_mm512_mul_pd(x1,x1);
    return 
    _mm512_mul_pd(x1,
        _mm512_fmadd_pd(x2,
            _mm512_fmadd_pd(x2,
                _mm512_fmadd_pd(x2,
                    _mm512_fmadd_pd(x2,
                        _mm512_fmadd_pd(x2,
                            _mm512_fmadd_pd(x2,
                                _mm512_fmadd_pd(x2,
                                    _mm512_fmadd_pd(x2,PCG32_AVX512_SET_DOUBLE_LLU(0x3CE880FF69A815D1LLU),PCG32_AVX512_SET_DOUBLE_LLU(0xBD6AE420DC08FB2ALLU)),
                                    PCG32_AVX512_SET_DOUBLE_LLU(0x3DE6123C686AD6A8LLU)
                                ),
                                PCG32_AVX512_SET_DOUBLE_LLU(0xBE5AE6454B5DC0B5LLU)
                            ),
                            PCG32_AVX512_SET_DOUBLE_LLU(0x3EC71DE3A524F063LL)
                        ),
                        PCG32_AVX512_SET_DOUBLE_LLU(0xBF2A01A01A013E1ALLU)
                    ),
                    PCG32_AVX512_SET_DOUBLE_LLU(0x3F811111111110B0LLU)
                ),
                PCG32_AVX512_SET_DOUBLE_LLU(0xBFC5555555555555LLU)
            ),
            PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU)
        )
    );
}

/**
 * 0.99999999999999998243004448007448662 + 
 * x2*(-0.499999999999999105881272803474436268 + 
 * x2*(0.041666666666658914344068844317924505 + 
 * x2*(-0.00138888888886231429175747130897185107 + 
 * x2*(0.0000248015872549765577961155967511699095 + 
 * x2*(-2.75573145508960795189972565635584642e-7 + 
 * x2*(2.08764776731016710219609723288490596e-9 + 
 * x2*(-1.14608862231521440480830153964369191e-11 + 
 * 4.58927688754481747776178904291483144e-14*x2)))))))
 */
static inline __m512d __x8__PCG32Cos(const __m512d x){
    const __m512d x2=_mm512_mul_pd(x,x);
    return 
    _mm512_fmadd_pd(x2,
        _mm512_fmadd_pd(x2,
            _mm512_fmadd_pd(x2,
                _mm512_fmadd_pd(x2,
                    _mm512_fmadd_pd(x2,
                        _mm512_fmadd_pd(x2,
                            _mm512_fmadd_pd(x2,
                                _mm512_fmadd_pd(x2,PCG32_AVX512_SET_DOUBLE_LLU(0x3D29D5D853164598LLU),PCG32_AVX512_SET_DOUBLE_LLU(0xBDA933E7C608646FLLU)),
                                PCG32_AVX512_SET_DOUBLE_LLU(0x3E21EEC9369F3E5DLL)
                            ),
                            PCG32_AVX512_SET_DOUBLE_LLU(0xBE927E4F82DB5BB1LLU)
                        ),
                        PCG32_AVX512_SET_DOUBLE_LLU(0x3EFA01A0192FB593LLU)
                    ),
                    PCG32_AVX512_SET_DOUBLE_LLU(0xBF56C16C16BF8D5DLLU)
                ),
                PCG32_AVX512_SET_DOUBLE_LLU(0x3FA55555555550F8LLU)
            ),
            PCG32_AVX512_SET_DOUBLE_LLU(0xBFDFFFFFFFFFFFF0LLU)
        ),
        PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU)
    );
}

// 0 < x <= 1
static inline __m512d __x8__PCG32Ln(const __m512d x){
    #if 1
    __m512d mantissa=_mm512_getmant_pd(x,_MM_MANT_NORM_1_2,_MM_MANT_SIGN_src);
    __m512d exponent=_mm512_getexp_pd(x);
    const __mmask8 mask=_mm512_cmp_pd_mask(mantissa,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF6A09E667F3BCDLLU),_CMP_GT_OQ);
    mantissa=_mm512_mask_mul_pd(mantissa,mask,mantissa,PCG32_AVX512_SET_DOUBLE_LLU(0x3FE0000000000000LLU));
    exponent=_mm512_mask_add_pd(exponent,mask,exponent,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU));
    const __m512d normalized=_mm512_sub_pd(mantissa,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU));
    const __m512d normalizedAdd2=_mm512_add_pd(mantissa,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU));
    const __m512d s=_mm512_div_pd(normalized,normalizedAdd2);
    const __m512d s2=_mm512_mul_pd(s,s);
    const __m512d s4=_mm512_mul_pd(s2,s2);
    static const long long unsigned int L[]={
        0x3FE5555555555593LLU,
        0x3FD999999997FA04LLU,
        0x3FD2492494229359LLU,
        0x3FCC71C51D8E78AFLLU,
        0x3FC7466496CB03DELLU,
        0x3FC39A09D078C69FLLU,
        0x3FC2F112DF3E5244LLU
    };
    const __m512d t1=_mm512_mul_pd(s2,_mm512_fmadd_pd(s4,_mm512_fmadd_pd(s4,_mm512_fmadd_pd(s4,PCG32_AVX512_SET_DOUBLE_LLU(L[6]),PCG32_AVX512_SET_DOUBLE_LLU(L[4])),PCG32_AVX512_SET_DOUBLE_LLU(L[2])),PCG32_AVX512_SET_DOUBLE_LLU(L[0])));
    const __m512d t2=_mm512_mul_pd(s4,_mm512_fmadd_pd(s4,_mm512_fmadd_pd(s4,PCG32_AVX512_SET_DOUBLE_LLU(L[5]),PCG32_AVX512_SET_DOUBLE_LLU(L[3])),PCG32_AVX512_SET_DOUBLE_LLU(L[1])));
    const __m512d R=_mm512_add_pd(t1,t2);
    const __m512d hfsq=_mm512_mul_pd(_mm512_mul_pd(normalized,normalized),PCG32_AVX512_SET_DOUBLE_LLU(0x3FE0000000000000LLU));
    const __m512d expln2High=_mm512_mul_pd(exponent,PCG32_AVX512_SET_DOUBLE_LLU(0x3FE62E42FEE00000LLU));
    const __m512d expln2Low =_mm512_mul_pd(exponent,PCG32_AVX512_SET_DOUBLE_LLU(0x3DEA39EF35793C76LLU));
    const __m512d result=_mm512_sub_pd(expln2High,_mm512_sub_pd(_mm512_sub_pd(hfsq,_mm512_fmadd_pd(s,_mm512_add_pd(hfsq,R),expln2Low)),normalized));
    #else
    const __m512i ix=_mm512_castpd_si512(x);
    const __m512d exponent=_mm512_getexp_pd(x);
    const __m512d mantissa=_mm512_getmant_pd(x,_MM_MANT_NORM_1_2,_MM_MANT_SIGN_src);
    const __m512d expLn2=_mm512_mul_pd(exponent,PCG32_AVX512_SET_DOUBLE_LLU(0x3FE62E42FEFA39EFLLU));
    const __m512d mantissaSub1=_mm512_sub_pd(mantissa,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU));
    const __m512d mantissaAdd1=_mm512_add_pd(mantissa,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU));
    const __m512d order=_mm512_div_pd(mantissaSub1,mantissaAdd1);
    const __m512d order2=_mm512_mul_pd(order,order);
    const long long unsigned int _2dive9=0x3FCC71C71C71C71CLLU;
    const long long unsigned int _2dive7=0x3FD2492492492492LLU;
    const long long unsigned int _2dive5=0x3FD999999999999ALLU;
    const long long unsigned int _2dive3=0x3FE5555555555555LLU;
    const long long unsigned int _2dive1=0x4000000000000000LLU;
    const __m512d result=_mm512_fmadd_pd(
        _mm512_fmadd_pd(
            _mm512_fmadd_pd(
                _mm512_fmadd_pd(
                    _mm512_fmadd_pd(
                        order2,PCG32_AVX512_SET_DOUBLE_LLU(_2dive9),PCG32_AVX512_SET_DOUBLE_LLU(_2dive7)),
                        order2,PCG32_AVX512_SET_DOUBLE_LLU(_2dive5)),
                        order2,PCG32_AVX512_SET_DOUBLE_LLU(_2dive3)),
                        order2,PCG32_AVX512_SET_DOUBLE_LLU(_2dive1)),
                        order,expLn2);
    #endif
    return result;
}

// 0 <= shift64 <= 7
static inline __m512d m512dRightShift(__m512d x,const unsigned shift64){
    static const __mmask8 mask[8]={
        (0b11111111<<0)&0xFF,
        (0b11111111<<1)&0xFF,
        (0b11111111<<2)&0xFF,
        (0b11111111<<3)&0xFF,
        (0b11111111<<4)&0xFF,
        (0b11111111<<5)&0xFF,
        (0b11111111<<6)&0xFF,
        (0b11111111<<7)&0xFF
    };
    #if 0
    static const __m512i m512dRightShiftIndexes[8]={
        _mm512_set_epi64(7,6,5,4,3,2,1,0),
        _mm512_set_epi64(6,5,4,3,2,1,0,0),
        _mm512_set_epi64(5,4,3,2,1,0,0,0),
        _mm512_set_epi64(4,3,2,1,0,0,0,0),
        _mm512_set_epi64(3,2,1,0,0,0,0,0),
        _mm512_set_epi64(2,1,0,0,0,0,0,0),
        _mm512_set_epi64(1,0,0,0,0,0,0,0),
        _mm512_set_epi64(0,0,0,0,0,0,0,0),
    };
    return _mm512_maskz_permutexvar_pd(mask[shift64],m512dRightShiftIndexes[shift64],x);
    #else
    switch(shift64){
    case 0:
        return x;
    case 1:
        return _mm512_maskz_permutexvar_pd(mask[1],_mm512_set_epi64(6,5,4,3,2,1,0,0),x);
    case 2:
        return _mm512_maskz_permutexvar_pd(mask[2],_mm512_set_epi64(5,4,3,2,1,0,0,0),x);
    case 3:
        return _mm512_maskz_permutexvar_pd(mask[3],_mm512_set_epi64(4,3,2,1,0,0,0,0),x);
    case 4:
        return _mm512_maskz_permutexvar_pd(mask[4],_mm512_set_epi64(3,2,1,0,0,0,0,0),x);
    case 5:
        return _mm512_maskz_permutexvar_pd(mask[5],_mm512_set_epi64(2,1,0,0,0,0,0,0),x);
    case 6:
        return _mm512_maskz_permutexvar_pd(mask[6],_mm512_set_epi64(1,0,0,0,0,0,0,0),x);
    case 7:
        return _mm512_maskz_permutexvar_pd(mask[7],_mm512_set_epi64(0,0,0,0,0,0,0,0),x);
    default:
        return x;
    }
    #endif
}

static inline void __x16__PCG32Normal(__x16__PCG32Struct* status,__x16__DoubleArray random,const double mu,const double sigma){
    #if 0

    __m512d u1=PCG32_AVX512_SET_DOUBLE_LLU(0x0000000000000000LLU);
    __m512d u2=PCG32_AVX512_SET_DOUBLE_LLU(0x0000000000000000LLU);
    __m512d S;
    unsigned validSaved=0;
    PCG32_AVX512_PCG32_INTRINSIC_LOAD_CONSTANT
    PCG32_AVX512_PCG32_INTRINSIC_CALCULATE
    PCG32_AVX512_UNIFORM_REAL_INTRINSIC
    PCG32_AVX512_UNIFORM_REAL_MIN_MAX_INTRINSIC(-1.0,1.0)
    const __m512d sum=_mm512_fmadd_pd(random1,random1,_mm512_mul_pd(random0,random0));
    const __mmask8 greaterThan0=_mm512_cmp_pd_mask(sum,PCG32_AVX512_SET_DOUBLE_LLU(0x0000000000000000LLU),_CMP_GT_OQ);
    const __mmask8 lessThan1   =_mm512_cmp_pd_mask(sum,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU),_CMP_LE_OQ);
    const __mmask8 valid=greaterThan0&lessThan1;
    unsigned validCount=__builtin_popcount(valid);
    if(validCount<8){
        validSaved=validCount;
        u1=_mm512_maskz_compress_pd(valid,random0);
        u2=_mm512_maskz_compress_pd(valid,random1);
        ReGenerate:
        PCG32_AVX512_PCG32_INTRINSIC_CALCULATE                                                              
        PCG32_AVX512_UNIFORM_REAL_INTRINSIC                                                
        PCG32_AVX512_UNIFORM_REAL_MIN_MAX_INTRINSIC(-1.0,1.0)
        const __m512d sum=_mm512_fmadd_pd(random1,random1,_mm512_mul_pd(random0,random0));
        const __mmask8 greaterThan0=_mm512_cmp_pd_mask(sum,PCG32_AVX512_SET_DOUBLE_LLU(0x0000000000000000LLU),_CMP_GT_OQ);
        const __mmask8 lessThan1   =_mm512_cmp_pd_mask(sum,PCG32_AVX512_SET_DOUBLE_LLU(0x3FF0000000000000LLU),_CMP_LE_OQ);
        const __mmask8 valid=greaterThan0&lessThan1;
        validCount=__builtin_popcount(valid);
        const __m512d compressed0=m512dRightShift(_mm512_maskz_compress_pd(valid,random0),validSaved);
        const __m512d compressed1=m512dRightShift(_mm512_maskz_compress_pd(valid,random1),validSaved);
        u1=_mm512_or_pd(u1,compressed0);
        u2=_mm512_or_pd(u2,compressed1);
        validSaved=validSaved+validCount;
        if(validSaved<8){
            goto ReGenerate;
        }
        S=_mm512_fmadd_pd(u1,u1,_mm512_mul_pd(u2,u2));
    }else{
        u1=random0;
        u2=random1;
        S=sum;
    }
    const __m512d lnS=__x8__PCG32Ln(S);
    const __m512d toMultiple=_mm512_sqrt_pd(_mm512_div_pd(_mm512_mul_pd(lnS,PCG32_AVX512_SET_DOUBLE_LLU(0xC000000000000000LLU)),S));
    __m512d result0=_mm512_mul_pd(toMultiple,u1);
    __m512d result1=_mm512_mul_pd(toMultiple,u2);
    
    #else

    PCG32_AVX512_PCG32_INTRINSIC_LOAD_CONSTANT
    __m512d negitive1=PCG32_AVX512_SET_DOUBLE_LLU(0xBFF0000000000000LLU);
    __m512d u1,u2;
    {
        PCG32_AVX512_PCG32_INTRINSIC_CALCULATE
        PCG32_AVX512_UNIFORM_REAL_INTRINSIC
        u1=_mm512_add_pd(random0,PCG32_AVX512_SET_DOUBLE_LLU(0x3DF0000000000000LLU));
        const long long unsigned int PiLLU      =0x400921FB54442D18LLU;
        const long long unsigned int _PiDive2LLU=0xBFF921FB54442D18LLU;
        random1=_mm512_mul_pd(random1,PCG32_AVX512_SET_DOUBLE_LLU(PiLLU));
        random1=_mm512_add_pd(random1,PCG32_AVX512_SET_DOUBLE_LLU(_PiDive2LLU));
        u2=random1;
    }
    __m512d Cos=__x8__PCG32Cos(u2);
    __m512d Sin=__x8__PCG32Sin(u2);
    {
        PCG32_AVX512_PCG32_INTRINSIC_CALCULATE
        __mmask8 mask=_mm512_test_epi64_mask(x,_mm512_set1_epi64(1));
        Cos=_mm512_mask_mul_pd(Cos,mask,Cos,negitive1);
        Sin=_mm512_mask_mul_pd(Sin,mask,Sin,negitive1);
    }
    __m512d toMultiple=_mm512_sqrt_pd(_mm512_mul_pd(__x8__PCG32Ln(u1),PCG32_AVX512_SET_DOUBLE_LLU(0xC000000000000000LLU)));
    __m512d result0=_mm512_mul_pd(toMultiple,Cos);
    __m512d result1=_mm512_mul_pd(toMultiple,Sin);

    #endif
    if(mu!=0.0||sigma!=1.0){
        union DoubleLLU{                                                                      
            double d;
            long long unsigned int u;
        }muUnion,sigmaUnion;
        muUnion.d   =mu;
        sigmaUnion.d=sigma;
        result0=_mm512_fmadd_pd(result0,PCG32_AVX512_SET_DOUBLE_LLU(sigmaUnion.u),PCG32_AVX512_SET_DOUBLE_LLU(muUnion.u));
        result1=_mm512_fmadd_pd(result1,PCG32_AVX512_SET_DOUBLE_LLU(sigmaUnion.u),PCG32_AVX512_SET_DOUBLE_LLU(muUnion.u));
    }
    _mm512_store_pd(random+0,result0);
    _mm512_store_pd(random+8,result1);
}

static inline void __x16__PCG32StandardNormal(__x16__PCG32Struct* status,__x16__DoubleArray random){
    __x16__PCG32Normal(status,random,0.0,1.0);
}

#endif

#ifdef __cplusplus
}
#endif

#if defined(__cplusplus)||PCG32_CUDA

template<typename Type>
PCG32_HOST_DEVICE static inline void PCG32SWAP(Type* array,const long long unsigned int index0,const long long unsigned int index1){
    const Type tempory=array[index0];                                                                            
    array[index0]=array[index1];                                                                                 
    array[index1]=tempory;  
}

template<typename Type>
PCG32_HOST_DEVICE static inline void PCG32UniformShuffle(PCG32Struct* status,Type* array,const long long unsigned int length){
    if(length>1){                                                                                
        for(long long unsigned int index=0;index<length-1;index=index+1){                        
            PCG32SWAP(array,index,PCG32Uniform_MaxBiggerThanMin(status,index,length-1));
        }                                                                                        
    }   
}

template<typename Type>
PCG32_HOST_DEVICE static inline void PCG32UniformShuffle_FirstK(PCG32Struct* status,Type* array,const long long unsigned int length,long long unsigned int k){
    if(length>1){                                                                                
        if(k>length){
            k=length;
        }
        for(long long unsigned int index=0;index<k;index=index+1){                        
            PCG32SWAP(array,index,PCG32Uniform_MaxBiggerThanMin(status,index,length-1));
        }                                                                                        
    }   
}

PCG32_HOST_DEVICE inline PCG32PRNG::result_type PCG32PRNG::operator()(){
    return Rand();
}

PCG32_HOST_DEVICE inline PCG32PRNG::PCG32PRNG():PCG32PRNG(0xADABF3924A46334BLLU){
}

PCG32_HOST_DEVICE inline PCG32PRNG::PCG32PRNG(long long unsigned int seed){
    SetSeed(seed);
}

PCG32_HOST_DEVICE inline void PCG32PRNG::SetSeed(long long unsigned int seed){
    PCG32SetSeed(&status,seed);
}

PCG32_HOST_DEVICE inline unsigned PCG32PRNG::Rand(){
    return PCG32(&status);
}

PCG32_HOST_DEVICE inline unsigned PCG32PRNG::Uniform(unsigned min,unsigned max){
    return PCG32Uniform(&status,min,max);
}

PCG32_HOST_DEVICE inline unsigned PCG32PRNG::Uniform_Strict(unsigned min,unsigned max){
    return PCG32Uniform_Strict(&status,min,max);
}

PCG32_HOST_DEVICE inline void PCG32PRNG::UniformSetStrictRange(const unsigned min,const unsigned max){
    PCG32UniformSetStrictRange(&status,min,max);
}

PCG32_HOST_DEVICE inline unsigned PCG32PRNG::Uniform_StrictRangeUnchanged(){
    return PCG32Uniform_StrictRangeUnchanged(&status);
}

PCG32_HOST_DEVICE inline unsigned PCG32PRNG::Uniform_MaxBiggerThanMin(const unsigned min,const unsigned max){
    return PCG32Uniform_MaxBiggerThanMin(&status,min,max);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::UniformReal(const double min,const double max){
    return PCG32UniformReal(&status,min,max);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::Normal(){
    return PCG32StandardNormal(&status);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::Normal(const double mu,const double sigma){
    return mu+Normal()*sigma;
}

PCG32_HOST_DEVICE inline bool PCG32PRNG::GammaInitialize(const double alpha,const double beta){
    return PCG32GammaInitialize(&status,alpha,beta);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::Gamma(){
    return PCG32Gamma(&status);
}

PCG32_HOST_DEVICE inline unsigned PCG32PRNG::Binomial(double probability,const unsigned repeatUnsigned){
    return PCG32Binomial(&status,probability,repeatUnsigned);
}

template<typename Type>
PCG32_HOST_DEVICE inline void PCG32PRNG::UniformShuffle(Type* array,const long long unsigned int length){
    PCG32UniformShuffle(&status,array,length);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::Exponential(const double lambda){
    return PCG32Exponential(&status,lambda);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::PowerLaw(const double min,const double alpha){
    return PCG32PowerLaw(&status,min,alpha);
}

PCG32_HOST_DEVICE inline void PCG32PRNG::RandomPointInCycle(const double radius,double* xy){
    PCG32RandomPointInCycle(&status,radius,xy);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::LogNormal(){
    return exp(Normal());
}

PCG32_HOST_DEVICE inline double PCG32PRNG::LogNormal(const double mu,const double sigma){
    return exp(Normal(mu,sigma));
}

PCG32_HOST_DEVICE inline unsigned PCG32PRNG::Geometric(const double probability){
    return PCG32Geometric(&status,probability);
}

PCG32_HOST_DEVICE inline bool PCG32PRNG::PoissonInitialize(const double mu){
    return PCG32PoissonInitialize(&status,mu);
}

PCG32_HOST_DEVICE inline double PCG32PRNG::Poisson(){
    return PCG32Poisson(&status);
}

#else

#define GENERATE_FOR_TYPE(TypeName,Type)                                                                         \
static inline void PCG32SWAP_##TypeName(Type* array,long long unsigned int index0,long long unsigned int index1){\
    const Type tempory=array[index0];                                                                            \
    array[index0]=array[index1];                                                                                 \
    array[index1]=tempory;                                                                                       \
}

GENERATE_FOR_TYPE(unsigned_char,unsigned char)
GENERATE_FOR_TYPE(unsigned_short,unsigned short)
GENERATE_FOR_TYPE(unsigned_long,unsigned long)
GENERATE_FOR_TYPE(char,char)
GENERATE_FOR_TYPE(short,short)
GENERATE_FOR_TYPE(int,int)
GENERATE_FOR_TYPE(long,long)
GENERATE_FOR_TYPE(float,float)
GENERATE_FOR_TYPE(double,double)
GENERATE_FOR_TYPE(unsigned_long_long,unsigned long long)
GENERATE_FOR_TYPE(long_long,long long)
GENERATE_FOR_TYPE(unsigned,unsigned)

#define PCG32_GENERIC_SWAP(array,index0,index1)           \
    _Generic((array),                                     \
        unsigned char*:      PCG32SWAP_unsigned_char,     \
        unsigned short*:     PCG32SWAP_unsigned_short,    \
        unsigned long*:      PCG32SWAP_unsigned_long,     \
        char*:               PCG32SWAP_char,              \
        short*:              PCG32SWAP_short,             \
        int*:                PCG32SWAP_int,               \
        long*:               PCG32SWAP_long,              \
        float*:              PCG32SWAP_float,             \
        double*:             PCG32SWAP_double,            \
        unsigned long long*: PCG32SWAP_unsigned_long_long,\
        long long*:          PCG32SWAP_long_long,         \
        default:             PCG32SWAP_unsigned           \
    )(array,index0,index1)

#define PCG32UniformShuffle(status,array,length)                                                     \
    do{                                                                                              \
        if(length>1){                                                                                \
            for(long long unsigned int index=0;index<length-1;index=index+1){                        \
                PCG32_GENERIC_SWAP(array,index,PCG32Uniform_MaxBiggerThanMin(status,index,length-1));\
            }                                                                                        \
        }                                                                                            \
    }while(0)

#define PCG32UniformShuffle_FirstK(status,array,length,k)                                            \
    do{                                                                                              \
        if(length>1){                                                                                \
            if(k>length){                                                                            \
                k=length;                                                                            \
            }                                                                                        \
            for(long long unsigned int index=0;index<k;index=index+1){                               \
                PCG32_GENERIC_SWAP(array,index,PCG32Uniform_MaxBiggerThanMin(status,index,length-1));\
            }                                                                                        \
        }                                                                                            \
    }while(0)

#endif

PCG32_HOST_DEVICE static inline unsigned rotr32(unsigned x,unsigned r){
    return x>>r|x<<(-r&31);
}

PCG32_HOST_DEVICE static inline long long unsigned int PCG32NextPrime(long long unsigned int prime){
    if(prime<=1){
        return 2;
    }
    if(prime==2){
        return 3;
    }
    if((prime&1)==0){
        prime=prime-1;
    }
    for(long long unsigned int number=prime+2;true;number=number+2){
        const long long unsigned int sqrtPrime=(long long unsigned int)(sqrt(number)+0.5);
        bool isPrime=true;
        for(long long unsigned int mod=3;mod<=sqrtPrime;mod=mod+2){
            if(number%mod==0){
                isPrime=false;
                break;
            }
        }
        if(isPrime){
            return number;
        }
    }
}

// This function requires that the inputs a and b be coprime.
PCG32_HOST_DEVICE static inline unsigned PCG32ModInverse(unsigned a,unsigned b){
    const long long m=b;
    long long x=1;
    long long lasxX=0;
    while(b){
        const long long quotient=a/b;
        long long temporary=a%b;
        a=b;
        b=temporary;
        temporary=x-quotient*lasxX;
        x=lasxX;
        lasxX=temporary;
    }
    if(a==1){
        return (x%m+m)%m;
    }
    return -1;
}

PCG32_HOST static inline long long unsigned int PCG32TimeNanoeconds(){
    struct timespec ts;
    timespec_get(&ts,TIME_UTC);
    return (long long unsigned int)ts.tv_sec*1000000000LLU+ts.tv_nsec;
}

// 0 < x
PCG32_HOST_DEVICE static inline double PCG32Ln(const double x){
    union{
        double d;
        long long unsigned int u;
    }u={x};
    const int exp2=((int)(u.u>>52)&0b011111111111)-1023;
    const double ln2=0.69314718055994530941723212145817656807550013436025525412068L;
    const double exp2ln2=exp2*ln2;
    u.u=u.u|(1023LLU<<52);
    const unsigned scale=(unsigned)((u.d-1)*20);
    static const double dive[20]={
                                                                   1.0L,
        0.952380952380952380952380952380952380952380952380952380952381L,
        0.909090909090909090909090909090909090909090909090909090909091L,
        0.869565217391304347826086956521739130434782608695652173913043L,
        0.833333333333333333333333333333333333333333333333333333333333L,
                                                                   0.8L,
        0.769230769230769230769230769230769230769230769230769230769231L,
        0.740740740740740740740740740740740740740740740740740740740741L,
        0.714285714285714285714285714285714285714285714285714285714286L,
        0.689655172413793103448275862068965517241379310344827586206897L,
        0.666666666666666666666666666666666666666666666666666666666667L,
        0.645161290322580645161290322580645161290322580645161290322581L,
                                                                 0.625L,
        0.606060606060606060606060606060606060606060606060606060606061L,
        0.588235294117647058823529411764705882352941176470588235294118L,
        0.571428571428571428571428571428571428571428571428571428571429L,
        0.555555555555555555555555555555555555555555555555555555555556L,
        0.540540540540540540540540540540540540540540540540540540540541L,
        0.526315789473684210526315789473684210526315789473684210526316L,
        0.512820512820512820512820512820512820512820512820512820512821L
    };
    static const double add[20]={
        0                                                               ,
        0.0487901641694320030653744042231646586079736644155824100400766L,
        0.0953101798043248600439521232807650922206053653086441991852398L,
         0.139761942375158697371529255667655342765778691851407511844627L,
         0.182321556793954626211718025154514633197389337914486983942726L,
         0.223143551314209755766295090309834503374601085548007213671288L,
         0.262364264467491052035495986880954397204166456131434140385718L,
         0.300104592450338080750512134625036338265870050479220125050075L,
         0.336472236621212930504593410216992090111483375313343466546742L,
         0.371563556432483033748048456219370829817911290933715848767662L,
         0.405465108108164381978013115464349136571990423462494197614014L,
          0.43825493093115525249394074839981643477333730749156374160271L,
         0.470003629245735553650937031148342064700899048812248040449392L,
         0.500775287912489242021965238745114228792595788771138396799254L,
          0.53062825106217039623154316318876232798710152395697181126391L,
          0.55961578793542268627088850052682659348608446086135068021803L,
         0.587786664902119008189731140618863769769379761376981181556741L,
         0.615185639090233450932872094888906388223475964178607933490905L,
         0.641853886172394775991035977203489329636277772670355842504632L,
         0.667829372575655434013509102345303533776156879593928337999732L
    };
    u.d=u.d*dive[scale];
    const double toAdd=exp2ln2+add[scale];
    const double order=(u.d-1.0)/(u.d+1.0);
    const double order2=order*order;
    return ((((2.0/9.0*order2+2.0/7.0)*order2+2.0/5.0)*order2+2.0/3.0)*order2+2.0)*order+toAdd;
}

PCG32_HOST_DEVICE static inline double PCG32PowerUnsigned(double x,unsigned n){
    double power=1.0;
    do{
        if(n&1){
            power=power*x;
        }
        n=n>>1;
        x=x*x;
    }while(n);
    return power;
}

PCG32_HOST_DEVICE static inline double PCG32LnStirlingLeft(const double x){
    static const double C[5]={1.0/12,-1.0/360,1.0/1260,-1.0/1680,1.0/1188};
    double xPower[5]; // x^1, x^3, x^5, x^7, x^9
    xPower[0]=1.0/x;
    const double x2=xPower[0]*xPower[0];
    xPower[1]=xPower[0]*x2;
    xPower[2]=xPower[1]*x2;
    xPower[3]=xPower[2]*x2;
    xPower[4]=xPower[3]*x2;
    return C[0]*xPower[0]+C[1]*xPower[1]+C[2]*xPower[2]+C[3]*xPower[3]+C[4]*xPower[4];
}

// only accept x > 0
PCG32_HOST_DEVICE static inline double PCG32LnStirling(double x){
    double sub=0.0;
    while(x<12){
        sub=sub+log(x);
        x=x+1;
    }
    const double ln2Pidive2=0.91893853320467274178032973640561763986139747363778341L;
    const double lnx=log(x);
    return ln2Pidive2+lnx*0.5+x*(lnx-1.0)+PCG32LnStirlingLeft(x)-sub;
}

#endif
