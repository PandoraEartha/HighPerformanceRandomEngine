# PCG32 - High-Performance Random Number Generator

A lightweight, high-performance pseudorandom number generator based on the PCG-XSH-RR algorithm, offering excellent statistical quality and computational efficiency. Provides both C and C++ interfaces with full CUDA compatibility for GPU acceleration, with optional AVX512 vectorization support.

# Performance

<img width="2144" height="2954" alt="pcg32_vs_std_speed_comparison_log" src="https://github.com/user-attachments/assets/18286bda-0793-419a-91b1-b8a16ce5465c" />

---

## 📊 Quick Reference

| Feature | Status |
|---------|--------|
| **Language** | C / C++ / CUDA |
| **License** | MIT |
| **Platforms** | Cross-platform (x86, x86_64, CUDA devices) |
| **AVX512** | Optional SIMD acceleration (16x parallel) |
| **Distributions** | 15+ statistical distributions |
| **Vectorization** | AVX512 batch processing (16 lanes) |

---

## ✨ Features

- ✅ **32-bit uniform random integers** `[0, 0xFFFFFFFF]`
- ✅ **Uniform real numbers** `[min, max)`
- ✅ **Multiple integer uniform variants** (general, strict, fast, preset-range)
- ✅ **Standard normal distribution** (Box-Muller transform)
- ✅ **Gamma distribution** (α ≥ 1, Marsaglia & Tsang method)
- ✅ **Binomial distribution** (efficient BTPE algorithm)
- ✅ **Poisson distribution** (saddlepoint approximation for large μ)
- ✅ **Exponential distribution**
- ✅ **Power-law (Pareto) distribution**
- ✅ **Geometric distribution** (standard + small-probability optimized)
- ✅ **Log-normal distribution** (standard + general)
- ✅ **Benford's law distributed numbers**
- ✅ **N-dimensional uniform points** in spheres
- ✅ **Uniform points** on/inside 2D circles
- ✅ **Uniform simplex sampling** (fixed-sum real variables)
- ✅ **Fisher-Yates array shuffling** (template for arbitrary types)
- ✅ **CUDA support** — all functions are `__host__ __device__`
- ✅ **AVX512 support** — generate 16 random numbers simultaneously
- ✅ **C++ class wrapper** — `PCG32PRNG` with convenient methods
- ✅ **C API** — prefix `PCG32` for all functions

---

## 🚀 Quick Start

### C++ Interface

```cpp
#include "PCG32.h"

// Basic initialization
PCG32PRNG rng(PCG32TimeNanoeconds(nullptr));

// Uniform real in [-1.0, 1.0)
double u = rng.UniformReal(-1.0, 1.0);

// Binomial distribution: 1000 trials with p=0.9
unsigned binom = rng.Binomial(0.9, 1000);

// Gamma distribution (must initialize first, alpha >= 1)
rng.GammaInitialize(2.0, 1.0);
double gamma = rng.Gamma();

// Fast strict-range integer generation
rng.UniformSetStrictRange(0, 999);
for (int i = 0; i < 100; ++i) {
    unsigned x = rng.Uniform_StrictRangeUnchanged();  // fast!
}

// Random integer (min/max auto-sorted)
unsigned x = rng.Uniform(10, 20);

// Shuffle an array
int myArray[100];
rng.UniformShuffle(myArray, 100);
```

### C Interface

```c
#include "PCG32.h"

PCG32Struct state;
PCG32SetSingleSeed(&state, PCG32TimeNanoeconds(NULL));

// Uniform real
double u = PCG32UniformReal(&state, -1.0, 1.0);

// Binomial
unsigned binom = PCG32Binomial(&state, 0.9, 1000);

// Gamma
PCG32GammaInitialize(&state, 2.0, 1.0);
double gamma = PCG32Gamma(&state);

// Fast strict-range preset
PCG32UniformSetStrictRange(&state, 0, 999);
unsigned x = PCG32Uniform_StrictRangeUnchanged(&state);

// Basic random integer
unsigned x = PCG32Uniform(&state, 10, 20);

// Shuffle array
int myArray[100];
PCG32UniformShuffle(&state, myArray, 100);
```

---

## 📚 API Reference

### Auxiliary Functions

| Function | Description |
|----------|-------------|
| `PCG32TimeNanoeconds` | Get current time in nanoseconds since epoch (1970-01-01 00:00:00 UTC) |

**Example:**
```c
long long unsigned int time = PCG32TimeNanoeconds();
PCG32Struct state;
PCG32SetSingleSeed(&state, time);
```

---

### Seed Initialization

| Function | Description |
|----------|-------------|
| `PCG32SetSingleSeed` | Initialize a single generator with a seed |
| `PCG32SetMultipleSeeds` | Initialize multiple generators with a seed array |
| `__x16__PCG32SetMultipleSeeds` | Initialize multiple AVX152 generators with a __x16__SeedArray array |

> ⚠️ **Note:** `PCG32SetSeed` is deprecated. Use `PCG32SetSingleSeed` for a single generator, or `PCG32SetMultipleSeeds` for multiple generators.

**Single Generator Example:**
```cpp
PCG32PRNG rng(PCG32TimeNanoeconds());
```

```c
PCG32Struct state;
PCG32SetSingleSeed(&state, PCG32TimeNanoeconds());
```

**Multi-threaded Seed Initialization Example:**
```cpp
#include "PCG32.h"
#include <omp.h>

#define THREAD_NUMBER 32

// or use malloc or std::vector if THREAD_NUMBER is not a compile-time constant
long long unsigned int seeds[THREAD_NUMBER]; 
PCG32Struct status[THREAD_NUMBER];
long long unsigned int time = PCG32TimeNanoeconds();

// Generate distinct seeds for each thread
for (unsigned threadIndex = 0; threadIndex < THREAD_NUMBER; threadIndex++) {
    seeds[threadIndex] = time + threadIndex * 0x123456789ABCDEFLLU + 1123;
}

// Initialize all generators with one call
PCG32SetMultipleSeeds(status, seeds, THREAD_NUMBER);

#pragma omp parallel num_threads(THREAD_NUMBER)
{
    int tid = omp_get_thread_num();
    // Each thread uses its own generator
    double u = PCG32UniformReal(&status[tid], 0.0, 1.0);
    // ... generate more numbers ...
}
```

**AVX512 Multi-threaded Seed Initialization Example:**

```c
#include "PCG32.h"
#include <omp.h>

#define THREAD_NUMBER 32

// or use malloc or std::vector if THREAD_NUMBER is not a compile-time constant
__x16__SeedArray AVX512Seeds[THREAD_NUMBER];
__x16__PCG32Struct AVX512Status[THREAD_NUMBER];
long long unsigned int time = PCG32TimeNanoeconds();

// Generate distinct seeds for each thread
for (unsigned threadIndex = 0; threadIndex < THREAD_NUMBER; threadIndex++) {
   for (unsigned index = 0; index < 16; index++) {
       AVX512Seeds[threadIndex][index] = time + threadIndex * 0x123456789ABCDEFLLU + index;
   }
}

// Initialize all generators for scalar and AVX512 usage
__x16__PCG32SetMultipleSeeds(AVX512Status, AVX512Seeds, THREAD_NUMBER);

#pragma omp parallel num_threads(THREAD_NUMBER)
{
   int tid = omp_get_thread_num();
   
   // AVX512 generator for this thread (generates 16 values at once)
   __x16__DoubleArray avxUniforms;
   __x16__PCG32UniformReal(&AVX512Status[tid], avxUniforms);
   // Process 16 random values...
}
```

**CUDA Seed Initialization Example:**
```cpp
#include "PCG32.h"

// CUDA kernel - each thread processes its own generator
__global__ void kernel(PCG32Struct* deviceStatus, unsigned* results, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    
    // Each thread copies its status from device memory to local stack
    // (register/local memory is much faster than global memory)
    PCG32Struct localStatus = deviceStatus[idx];
    
    // Generate random numbers using local state
    for (int i = 0; i < 10; i++) {
        unsigned random = PCG32(&localStatus);
        results[idx * 10 + i] = random;
    }
    
    // Write back updated state if needed for subsequent kernel launches
    deviceStatus[idx] = localStatus;
}

int main() {
    const int N = 1024 * 1024;  // 1 million threads
    const int BLOCK_SIZE = 256;
    const int GRID_SIZE = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    // Host-side arrays
    long long unsigned int* seeds = (long long unsigned int*)malloc(N * sizeof(long long unsigned int));
    PCG32Struct* hostStatus = (PCG32Struct*)malloc(N * sizeof(PCG32Struct));
    unsigned* hostResults = (unsigned*)malloc(N * 10 * sizeof(unsigned));
    
    // Device-side arrays
    PCG32Struct* deviceStatus;
    unsigned* deviceResults;
    cudaMalloc(&deviceStatus, N * sizeof(PCG32Struct));
    cudaMalloc(&deviceResults, N * 10 * sizeof(unsigned));
    
    // Generate seeds on host
    long long unsigned int time = PCG32TimeNanoeconds();
    for (int i = 0; i < N; i++) {
        seeds[i] = time + i * 0x123456789ABCDEFLLU + 1123;
    }
    
    // Initialize all generators on host using PCG32SetMultipleSeeds
    PCG32SetMultipleSeeds(hostStatus, seeds, N);
    
    // Copy initialized generators to device
    cudaMemcpy(deviceStatus, hostStatus, N * sizeof(PCG32Struct), cudaMemcpyHostToDevice);
    
    // Launch kernel with proper grid configuration
    kernel<<<GRID_SIZE, BLOCK_SIZE>>>(deviceStatus, deviceResults, N);
    cudaDeviceSynchronize();
    
    // Copy results back
    cudaMemcpy(hostResults, deviceResults, N * 10 * sizeof(unsigned), cudaMemcpyDeviceToHost);
    
    // Clean up
    free(seeds);
    free(hostStatus);
    free(hostResults);
    cudaFree(deviceStatus);
    cudaFree(deviceResults);
    
    return 0;
}
```

---

### Core Functions

| Function | Description |
|----------|-------------|
| `PCG32` | Generate 32-bit uniform integer `[0, 0xFFFFFFFF]` |

---

### Uniform Integer Distributions

| Function | Description |
|----------|-------------|
| `PCG32Uniform` | Uniform integer `[min, max]` (auto-orders min/max) |
| `PCG32Uniform_Strict` | Uniform integer **requires** `min <= max` and gap not power-of-2 |
| `PCG32UniformSetStrictRange` | Preset strict range for repeated use |
| `PCG32Uniform_StrictRangeUnchanged` | Generate using preset strict range (fastest!) |
| `PCG32Uniform_MaxBiggerThanMin` | Fast uniform (assumes `min <= max`) |

**Example:**
```cpp
// Auto-order (safe)
unsigned x1 = rng.Uniform(20, 10);  // min=10, max=20

// Strict (faster, but you guarantee min <= max)
unsigned x2 = rng.Uniform_Strict(0, 999);

// Preset range (fastest for repeated use)
rng.UniformSetStrictRange(0, 999);
for (int i = 0; i < 1000; ++i) {
    unsigned x3 = rng.Uniform_StrictRangeUnchanged();  // no range checks!
}
```

```c
// C API equivalent
PCG32UniformSetStrictRange(&state, 0, 999);
unsigned x3 = PCG32Uniform_StrictRangeUnchanged(&state);
```

---

### Uniform Real Distributions

| Function | Description |
|----------|-------------|
| `PCG32UniformReal` | Uniform real in `[min, max)` |

**Example:**
```cpp
double u1 = rng.UniformReal(0.0, 1.0);     // [0, 1)
double u2 = rng.UniformReal(-10.0, 10.0);  // [-10, 10)
```

```c
double u = PCG32UniformReal(&state, 0.0, 1.0);
```

---

### Normal Distribution

| Function | Description |
|----------|-------------|
| `PCG32StandardNormal` | Standard normal `N(0,1)` |
| `PCG32StandardNormal2D` | Generate 2 standard normal samples at once |
| `PCG32StandardNormal3D` | Generate 3 standard normal samples at once |
| `PCG32StandardNormalNDimension` | Generate N standard normal samples |

**Example:**
```cpp
double z = rng.Normal();              // N(0,1)
double x = rng.Normal(5.0, 2.0);      // N(5, 4)
```

```c
double z = PCG32StandardNormal(&state);
```

---

### Gamma Distribution

| Function | Description |
|----------|-------------|
| `PCG32GammaInitialize` | Initialize Gamma distribution (α ≥ 1) |
| `PCG32Gamma` | Generate `Gamma(α, β)` sample |

**Example:**
```cpp
if (rng.GammaInitialize(2.0, 1.0)) {  // α=2.0, β=1.0
    double gamma = rng.Gamma();
}
```

```c
if (PCG32GammaInitialize(&state, 2.0, 1.0)) {
    double gamma = PCG32Gamma(&state);
}
```

> ⚠️ **Note:** Only supports shape parameter `α >= 1`.

---

### Binomial Distribution

| Function | Description |
|----------|-------------|
| `PCG32Binomial` | Binomial(n, p) |

**Example:**
```cpp
// 1000 trials, success probability 0.9
unsigned successes = rng.Binomial(0.9, 1000);
```

```c
unsigned successes = PCG32Binomial(&state, 0.9, 1000);
```

---

### Poisson Distribution

| Function | Description |
|----------|-------------|
| `PCG32PoissonInitialize` | Initialize Poisson distribution |
| `PCG32Poisson` | Generate Poisson(μ) sample |

**Example:**
```cpp
if (rng.PoissonInitialize(5.0)) {
    double poisson = rng.Poisson();
}
```

```c
if (PCG32PoissonInitialize(&state, 5.0)) {
    double poisson = PCG32Poisson(&state);
}
```

---

### Exponential Distribution

| Function | Description |
|----------|-------------|
| `PCG32Exponential` | Exponential(λ) |

**Example:**
```cpp
double exp = rng.Exponential(2.0);  // λ = 2.0
```

```c
double exp = PCG32Exponential(&state, 2.0);
```

---

### Power-Law Distribution

| Function | Description |
|----------|-------------|
| `PCG32PowerLaw` | Power-law (min, α) |

**Example:**
```cpp
double pl = rng.PowerLaw(1.0, 2.5);  // min=1.0, α=2.5
```

```c
double pl = PCG32PowerLaw(&state, 1.0, 2.5);
```

---

### Geometric Distribution

| Function | Description |
|----------|-------------|
| `PCG32Geometric` | Geometric(p) |
| `PCG32Geometric_SmallProbability` | Geometric(p) optimized for small p |

**Example:**
```cpp
unsigned g1 = rng.Geometric(0.5);                 // standard
unsigned g2 = rng.Geometric(0.001);               // uses standard
// For very small p, use specialized version:
unsigned g3 = PCG32Geometric_SmallProbability(&state, 0.0001);
```

---

### Log-Normal Distribution

| Function | Description |
|----------|-------------|
| `PCG32StandardLogNormal` | Standard log-normal (μ=0, σ=1) |
| `PCG32LogNormal` | Log-normal(μ, σ) |

**Example:**
```cpp
double ln1 = rng.LogNormal();              // μ=0, σ=1
double ln2 = rng.LogNormal(1.0, 0.5);      // μ=1, σ=0.5
```

```c
double ln1 = PCG32StandardLogNormal(&state);
double ln2 = PCG32LogNormal(&state, 1.0, 0.5);
```

---

### Benford's Law

| Function | Description |
|----------|-------------|
| `PCG32Benford` | Benford's law (digits 1-12) |
| `PCG32Benford_SpecificLength` | Benford with specified digit length range |

**Example:**
```cpp
double benford = rng.Benford();                     // 1-12 digits
double benford2 = rng.Benford_SpecificLength(3, 5); // 3-5 digits
```

```c
double benford = PCG32Benford(&state);
double benford2 = PCG32Benford_SpecificLength(&state, 3, 5);
```

---

### Geometric Functions

#### 3D Sphere

| Function | Description |
|----------|-------------|
| `PCG32RandomPointInSphere3D` | Uniform point in 3D sphere |

**Example:**
```cpp
double point[3];
rng.RandomPointInSphere3D(5.0, point);  // radius=5.0
// point[0], point[1], point[2] are uniformly distributed
```

```c
double point[3];
PCG32RandomPointInSphere3D(&state, 5.0, point);
```

#### N-Dimensional Sphere

| Function | Description |
|----------|-------------|
| `PCG32RandomPointInSphereNDimension` | Uniform point in N-dimensional sphere |

**Example:**
```cpp
double point[10];
rng.RandomPointInSphereNDimension(3.0, 10, point);  // radius=3.0, dim=10
```

```c
double point[10];
PCG32RandomPointInSphereNDimension(&state, 3.0, 10, point);
```

#### 2D Circle

| Function | Description |
|----------|-------------|
| `PCG32RandomPointInCycle` | Uniform point in 2D circle |

**Example:**
```cpp
double xy[2];
rng.RandomPointInCycle(2.0, xy);  // radius=2.0
// xy[0], xy[1] uniformly distributed inside circle
```

```c
double xy[2];
PCG32RandomPointInCycle(&state, 2.0, xy);
```

---

### Simplex Sampling

| Function | Description |
|----------|-------------|
| `PCG32UniformSumReal` | N uniform reals summing to a fixed value |

**Example:**
```cpp
double vars[5];
rng.UniformSumReal(5, 10.0, vars);  // 5 variables summing to 10.0
```

```c
double vars[5];
PCG32UniformSumReal(&state, 5, 10.0, vars);
```

---

### Shuffling

| Function | Description |
|----------|-------------|
| `PCG32UniformShuffle` | Fisher-Yates shuffle (type-generic) |
| `PCG32UniformShuffle_FirstK` | Fisher-Yates shuffle first K elements |

**C++ Example:**
```cpp
int array[100];
for (int i = 0; i < 100; ++i) array[i] = i;
rng.UniformShuffle(array, 100);              // shuffle all
rng.UniformShuffle_FirstK(array, 100, 20);   // shuffle first 20 only
```

**C Example:**
```c
int array[100];
for (int i = 0; i < 100; ++i) array[i] = i;
PCG32UniformShuffle(&state, array, 100);
PCG32UniformShuffle_FirstK(&state, array, 100, 20);
```

> ⚠️ **C API Note:** The C macro supports basic types (int, float, double, etc.). For custom types, use the C++ interface.

---

## ⚡ AVX512 Vectorization

When compiled with GCC/G++ on x86_64 with AVX512F and AVX512DQ, the library provides SIMD functions that generate **16 random numbers simultaneously**.

### Build Example

```bash
g++ -O3 -mavx512f -mavx512dq -mfma -std=c++11 -o myapp main.cpp
```

### Data Types (64-byte aligned)

| Type | Description |
|------|-------------|
| `__x16__StateArray` | 16 × uint64_t (generator states) |
| `__x16__SeedArray` | 16 × uint64_t (seeds) |
| `__x16__UnsignedArray` | 16 × uint32_t (random integers) |
| `__x16__DoubleArray` | 16 × double (random reals) |

### AVX512 Functions

| Function | Description |
|----------|-------------|
| `__x16__PCG32SetSingleSeed` | Initialize 16 generators with a seed array |
| `__x16__PCG32SetMultipleSeeds` | Initialize multiple 16-lane generator sets |
| `__x16__PCG32` | Generate 16 random 32-bit integers |
| `__x16__PCG32UniformReal` | Generate 16 uniform reals in `[0,1)` |
| `__x16__PCG32UniformReal_MinMax` | Generate 16 uniform reals in `[min, max)` |
| `__x16__PCG32UniformSetStrictRange` | Preset a strict range for all 16 lanes |
| `__x16__PCG32Uniform_StrictRangeUnchanged` | Generate 16 integers from preset range |
| `__x16__PCG32StandardNormal` | Generate 16 standard normal N(0,1) samples |
| `__x16__PCG32Normal` | Generate 16 normal N(mu, sigma) samples |

> ⚠️ **Note:** `__x16__PCG32SetSeed` is no longer support. Use `__x16__PCG32SetSingleSeed` for a single 16-lane set, or `__x16__PCG32SetMultipleSeeds` for multiple sets.

---

### AVX512 Example

```cpp
#include "PCG32.h"

// 1. __x16__PCG32SetSingleSeed - Initialize 16 generators with seeds
__x16__PCG32Struct avxState;
__x16__SeedArray seeds = {
    0x123456789ABCDEF0ULL, 0x23456789ABCDEF01ULL,  // ... 16 seeds total
    // ... fill all 16 seeds
};
__x16__PCG32SetSingleSeed(&avxState, seeds);

// 2. __x16__PCG32 - Generate 16 random 32-bit integers
__x16__UnsignedArray randoms;
__x16__PCG32(&avxState, randoms);
// randoms[0] through randoms[15] contain random 32-bit values

// 3. __x16__PCG32UniformReal - Generate 16 uniform reals in [0, 1)
__x16__DoubleArray uniforms;
__x16__PCG32UniformReal(&avxState, uniforms);
// uniforms[0] through uniforms[15] are in [0.0, 1.0)

// 4. __x16__PCG32UniformSetStrictRange - Preset strict range for all lanes
//    Range: [0, 999] for all 16 generators
__x16__PCG32UniformSetStrictRange(&avxState, 0, 999);

// 5. __x16__PCG32Uniform_StrictRangeUnchanged - Generate 16 integers from preset range
__x16__UnsignedArray strictRandoms;
for (int batch = 0; batch < 10; ++batch) {
    __x16__PCG32Uniform_StrictRangeUnchanged(&avxState, strictRandoms);
    // strictRandoms[0..15] are all in [0, 999]
    // Process the 16 random values...
}
// Note: For optimal performance with AVX512 strict range, the range should be
// preset once with __x16__PCG32UniformSetStrictRange, then reuse it with
// __x16__PCG32Uniform_StrictRangeUnchanged for many batches.

// 6. Generate 16 uniform reals in custom range [min, max)
__x16__DoubleArray customUniforms;
__x16__PCG32UniformReal_MinMax(&avxState, -5.0, 5.0, customUniforms);
// customUniforms[0..15] are all in [-5.0, 5.0)

// 7. Generate 16 standard normal N(0,1) samples
__x16__DoubleArray normalSamples;
__x16__PCG32StandardNormal(&avxState, normalSamples);
// normalSamples[0..15] are standard normal distributed

// 8. Generate 16 normal N(mu, sigma) samples
__x16__DoubleArray customNormalSamples;
__x16__PCG32Normal(&avxState, customNormalSamples, 2.0, 0.5);
// customNormalSamples[0..15] are N(2.0, 0.5) distributed
```

---

### AVX512 Multi-threaded Seed Initialization

```cpp
#include "PCG32.h"
#include <omp.h>

#define THREAD_NUMBER 32

// or use malloc or std::vector if THREAD_NUMBER is not a compile-time constant
__x16__SeedArray AVX512Seeds[THREAD_NUMBER];
__x16__PCG32Struct AVX512Status[THREAD_NUMBER];
long long unsigned int time = PCG32TimeNanoeconds();

// Generate distinct seeds for each thread
for (unsigned threadIndex = 0; threadIndex < THREAD_NUMBER; threadIndex++) {
    for (unsigned index = 0; index < 16; index++) {
        AVX512Seeds[threadIndex][index] = time + threadIndex * 0x123456789ABCDEFLLU + index;
    }
}

// Initialize all generators for scalar and AVX512 usage
__x16__PCG32SetMultipleSeeds(AVX512Status, AVX512Seeds, THREAD_NUMBER);

#pragma omp parallel num_threads(THREAD_NUMBER)
{
    int tid = omp_get_thread_num();
    
    // AVX512 generator for this thread (generates 16 values at once)
    __x16__DoubleArray avxUniforms;
    __x16__PCG32UniformReal(&AVX512Status[tid], avxUniforms);
    // Process 16 random values...
}
```

> ⚠️ **Note:** AVX512 functions are not available in CUDA mode.

---

## 🎮 CUDA Support

All functions are decorated with `PCG32_HOST_DEVICE`, enabling direct usage within CUDA kernel code. Simply include this header in your `.cu` files.

### CUDA Example

```cpp
#include "PCG32.h"

// CUDA kernel - each thread processes its own generator
__global__ void kernel(PCG32Struct* deviceStatus, unsigned* results, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    
    // Each thread copies its status from device memory to local stack
    // (register/local memory is much faster than global memory)
    PCG32Struct localStatus = deviceStatus[idx];
    
    // Generate random numbers using local state
    for (int i = 0; i < 10; i++) {
        unsigned random = PCG32(&localStatus);
        results[idx * 10 + i] = random;
    }
    
    // Write back updated state if needed for subsequent kernel launches
    deviceStatus[idx] = localStatus;
}

int main() {
    const int N = 1024 * 1024;  // 1 million threads
    const int BLOCK_SIZE = 256;
    const int GRID_SIZE = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    // Host-side arrays
    long long unsigned int* seeds = (long long unsigned int*)malloc(N * sizeof(long long unsigned int));
    PCG32Struct* hostStatus = (PCG32Struct*)malloc(N * sizeof(PCG32Struct));
    unsigned* hostResults = (unsigned*)malloc(N * 10 * sizeof(unsigned));
    
    // Device-side arrays
    PCG32Struct* deviceStatus;
    unsigned* deviceResults;
    cudaMalloc(&deviceStatus, N * sizeof(PCG32Struct));
    cudaMalloc(&deviceResults, N * 10 * sizeof(unsigned));
    
    // Generate seeds on host
    long long unsigned int time = PCG32TimeNanoeconds();
    for (int i = 0; i < N; i++) {
        seeds[i] = time + i * 0x123456789ABCDEFLLU + 1123;
    }
    
    // Initialize all generators on host using PCG32SetMultipleSeeds
    PCG32SetMultipleSeeds(hostStatus, seeds, N);
    
    // Copy initialized generators to device
    cudaMemcpy(deviceStatus, hostStatus, N * sizeof(PCG32Struct), cudaMemcpyHostToDevice);
    
    // Launch kernel with proper grid configuration
    kernel<<<GRID_SIZE, BLOCK_SIZE>>>(deviceStatus, deviceResults, N);
    cudaDeviceSynchronize();
    
    // Copy results back
    cudaMemcpy(hostResults, deviceResults, N * 10 * sizeof(unsigned), cudaMemcpyDeviceToHost);
    
    // Clean up
    free(seeds);
    free(hostStatus);
    free(hostResults);
    cudaFree(deviceStatus);
    cudaFree(deviceResults);
    
    return 0;
}
```

---

## 🛠️ Build Instructions

### Standard Build

```bash
# C++
g++ -O3 -std=c++11 -o myapp main.cpp

# C
gcc -O3 -std=c11 -o myapp main.c
```

### CUDA Build

```bash
nvcc -O3 -std=c++11 -arch=sm_70 -o myapp main.cu
```

### AVX512 Build

```bash
# GCC/G++ only
g++ -O3 -mavx512f -mavx512dq -mfma -std=c++11 -o myapp main.cpp
```

---

## 📋 Complete Distribution Quick Reference

| Distribution | C Function | C++ Method |
|--------------|------------|------------|
| Uniform Integer | `PCG32Uniform` | `Uniform` |
| Uniform Strict | `PCG32Uniform_Strict` | `Uniform_Strict` |
| Uniform Preset | `PCG32Uniform_StrictRangeUnchanged` | `Uniform_StrictRangeUnchanged` |
| Uniform Real | `PCG32UniformReal` | `UniformReal` |
| Normal | `PCG32StandardNormal` | `Normal` |
| Gamma | `PCG32Gamma` | `Gamma` |
| Binomial | `PCG32Binomial` | `Binomial` |
| Poisson | `PCG32Poisson` | `Poisson` |
| Exponential | `PCG32Exponential` | `Exponential` |
| Power-Law | `PCG32PowerLaw` | `PowerLaw` |
| Geometric | `PCG32Geometric` | `Geometric` |
| Log-Normal | `PCG32LogNormal` | `LogNormal` |
| Benford | `PCG32Benford` | `Benford` |

---

## ⚠️ Important Notes

1. **Always set a seed** before generating numbers. Use `PCG32SetSingleSeed` for a single generator, or `PCG32SetMultipleSeeds` for multiple generators. `PCG32SetSeed` is deprecated.
2. **Multi-threading**: Use different `PCG32Struct` per thread with different seeds.
3. **Gamma restriction**: Currently only supports shape parameter `α >= 1`.
4. **Performance tip**: For repeated generation within the same integer range, use `UniformSetStrictRange` + `Uniform_StrictRangeUnchanged` for optimal speed.
5. **AVX512 restriction**: Only available with GCC/G++ on x86_64 with `-mavx512f -mavx512dq -mfma` flags. Not available in CUDA mode.
6. **CUDA**: Set seeds at host side using PCG32SetMultipleSeeds. All functions are `__host__ __device__` compatible.

---

## 📖 Documentation

For complete documentation, see the header file comments or visit:
- [GitHub Repository](https://github.com/PandoraEartha/HighPerformanceRandomEngine)

---

## 📝 License

MIT License — feel free to use in commercial and open-source projects.

---

## 👤 Author

**PandoraEartha**

- GitHub: [@PandoraEartha](https://github.com/PandoraEartha)
- Project: [HighPerformanceRandomEngine](https://github.com/PandoraEartha/HighPerformanceRandomEngine)

---

## 🙏 Acknowledgments

- PCG algorithm by M.E. O'Neill
- Marsaglia & Tsang Gamma method
- BTPE algorithm for Binomial distribution
