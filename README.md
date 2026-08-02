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
PCG32PRNG rng(time(nullptr));

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
PCG32SetSeed(&state, time(NULL));

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

### Core Functions

| Function | Description |
|----------|-------------|
| `PCG32SetSeed` | Initialize generator with seed |
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

### Build Flags

```bash
g++ -O3 -mavx512f -mavx512dq -std=c++11 -o myapp main.cpp
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
| `__x16__PCG32SetSeed` | Initialize 16 generators |
| `__x16__PCG32` | Generate 16 random 32-bit integers |
| `__x16__PCG32UniformReal` | Generate 16 uniform reals in `[0,1)` |
| `__x16__PCG32UniformSetStrictRange` | Preset a strict range for all 16 lanes |
| `__x16__PCG32Uniform_StrictRangeUnchanged` | Generate 16 integers from preset range |

### Example

```cpp
#include "PCG32.h"

// 1. Initialize 16 generators with seeds
__x16__PCG32Struct avxState;
__x16__SeedArray seeds = {
    0x123456789ABCDEF0ULL, 0x23456789ABCDEF01ULL,
    0x3456789ABCDEF012ULL, 0x456789ABCDEF0123ULL,
    0x56789ABCDEF01234ULL, 0x6789ABCDEF012345ULL,
    0x789ABCDEF0123456ULL, 0x89ABCDEF01234567ULL,
    0x9ABCDEF012345678ULL, 0xABCDEF0123456789ULL,
    0xBCDEF0123456789AULL, 0xCDEF0123456789ABULL,
    0xDEF0123456789ABCULL, 0xEF0123456789ABCDULL,
    0xF0123456789ABCDEULL, 0x0123456789ABCDEFULL
};
__x16__PCG32SetSeed(&avxState, seeds);

// 2. Generate 16 random integers
__x16__UnsignedArray randoms;
__x16__PCG32(&avxState, randoms);
// randoms[0] through randoms[15] contain random 32-bit values

// 3. Generate 16 uniform reals in [0, 1)
__x16__DoubleArray uniforms;
__x16__PCG32UniformReal(&avxState, uniforms);

// 4. Preset strict range for all lanes
__x16__PCG32UniformSetStrictRange(&avxState, 0, 999);

// 5. Generate 16 integers from preset range (repeated use for performance)
for (int batch = 0; batch < 10; ++batch) {
    __x16__UnsignedArray strictRandoms;
    __x16__PCG32Uniform_StrictRangeUnchanged(&avxState, strictRandoms);
    // strictRandoms[0..15] are all in [0, 999]
}
```

> ⚠️ **Note:** AVX512 functions are not available in CUDA mode.

---

## 🎮 CUDA Support

All functions are decorated with `PCG32_HOST_DEVICE`, enabling direct usage within CUDA kernel code. Simply include this header in your `.cu` files.

### CUDA Example

```cpp
#include "PCG32.h"

__global__ void kernel(PCG32Struct* states, int N, double* results) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    
    PCG32Struct* state = &states[idx];
    double u = PCG32UniformReal(state, 0.0, 1.0);
    results[idx] = u;
}

int main() {
    const int N = 1000;
    PCG32Struct* d_states;
    cudaMalloc(&d_states, N * sizeof(PCG32Struct));
    
    // Initialize each state with a unique seed
    // ... (initialize on device)
    
    kernel<<<N/256, 256>>>(d_states, N, d_results);
    // ...
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
g++ -O3 -mavx512f -mavx512dq -std=c++11 -o myapp main.cpp
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

1. **Always set a seed** before generating numbers.
2. **Multi-threading**: Use separate generator instances per thread with distinct seeds.
3. **Gamma restriction**: Currently only supports shape parameter `α >= 1`.
4. **Performance tip**: For repeated generation within the same integer range, use `UniformSetStrictRange` + `Uniform_StrictRangeUnchanged` for optimal speed.
5. **AVX512 restriction**: Only available with GCC/G++ on x86_64 with `-mavx512f -mavx512dq` flags. Not available in CUDA mode.
6. **CUDA**: All functions are `__host__ __device__` compatible.

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
