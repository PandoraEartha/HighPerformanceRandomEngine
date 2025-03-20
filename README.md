[中文介绍](https://zhuanlan.zhihu.com/p/27708881274)

# High Performance C Random Engine
Head only high performance pseudorandom engine base on PCG-XSH-RR. Multithreading allowed, easy to use, statistically good. 

仅有头文件的随机引擎, 基于PCG-XSH-RR算法, 允许多线程, 简单易用, 随机性强. 

## How to use

`#include "PCG32.h"`

```
PCG32Struct PCGStatus;
PCG32SetSeed(&PCGStatus,time(NULL));
double random=PCG32UniformReal(&PCGStatus,-1,999);
```

note:

1.Remember to set the seed. 记得设置种子. 

2.Use its own `PCG32Struct` in each thread function and set different seed. 每个线程函数采用各自的`PCG32Struct`并且设置单独的种子

# Functions 

`void PCG32SetSeed(PCG32Struct* status,long long unsigned int seed);`

Set the seed of PCG random engine and initlize 

初始化随机引擎, 设置随机种子

`unsigned PCG32(PCG32Struct* status);`

Generate a unsigned type random number in range of [0,0xFFFFFFFF(4294967295)]

产生[0,0xFFFFFFFF(4294967295)]范围的`unsigned`类型的随机整数

`unsigned PCG32Uniform(PCG32Struct* status,unsigned min,unsigned max);`

Generate a unsigned type random number that obey uniform distrubution in range of [min,max]

产生[min,max]范围的`unsigned`类型的随机整数

`double PCG32UniformReal(PCG32Struct* status,double min,double max);`

Generate a double type random number that obey uniform distrubution in range of [min,max]

产生[min,max]范围的`double`类型的随机整数

`double PCG32StandardNormal(PCG32Struct* status);`

Generate a double type random number that obey standard normal distrubution

产生`double`类型的服从标准正态分布的随机数

`void PCG32UniformShuffle(PCG32Struct* status,unsigned* array,unsigned length);`

Shuffling, swaps elements in the input sequence to ensure that every element etc. can appear anywhere in the array.

洗牌, 对输入序列进行元素交换, 保证每个元素等可能的出现在数组中任意位置. 

`unsigned PCG32Uniform_Strict(PCG32Struct* status,const unsigned min,const unsigned max);`

Generate a unsigned type random number that obey uniform distrubution in range of [min,max], `max` can not smaller than `min` and gap(`max`-`min`) can not be power of 2

产生[min,max]范围的`unsigned`类型的随机整数, `max`不能小于`min`且gap(`max`-`min`)不能是2的幂


# Performance

Performance test base on 13490F WSL2 Ubuntu, built using g++ 14.02. Compare with `std::default_random_engine` in C++ `random`. The result below shows the time consuming on generating 4074074037 random number on spicific distribution.

性能测试基于13490F WSL2 Ubuntu, 采用g++ 14.02编译. 与C++`random`库中的`std::default_random_engine`进行比较. 下面的结果显示了在特定分布上生成4074074037随机数所花费的时间. 

![image](https://github.com/user-attachments/assets/80df23e7-8cc0-43ed-885a-6817652c5654)

## Uniform Real Distribution 

![image](https://github.com/user-attachments/assets/37c8ee99-9fec-4b91-8cc9-c41bf580cf78)

`PCG32UniformReal` is 8.6399 times faster than `std::uniform_real_distribution`

`PCG32UniformReal` 比 `std::uniform_real_distribution` 快8.6399倍

## Normal Distribution 

![image](https://github.com/user-attachments/assets/a0038c41-7b70-4ec0-a77e-4215777b70dc)

`PCG32StandardNormal` is 5.8596 times faster than `std::normal_distribution`

`PCG32StandardNormal` 比 `std::normal_distribution` 快5.8596倍

## Base Generator

![image_2025-03-04_21-08-13](https://github.com/user-attachments/assets/66eb203c-1693-461d-9fa4-733996eb8351)

`PCG32` is 1.9767 times faster than `std::default_random_engine`

`PCG32` 比`std::default_random_engine`快1.9898倍

## Uniform Unsigned

![image_2025-03-04_21-10-18](https://github.com/user-attachments/assets/30c4da1d-7ad5-47da-a976-3c3783f3ac56)

`PCG32Uniform` is 1.7713 times faster than `std::uniform_int_distribution<unsigned>`

`PCG32Uniform` 比 `std::uniform_int_distribution<unsigned>` 快1.7713倍

# Test Demo

`g++ -o test test.cpp -O2 && ./test`

## `PCG32Test.h`

## Functions 

`double KolmogorovSmirnovTest(double* data,long long unsigned int length,double (*CDF)(double));`

`double KSTest(double* data,long long unsigned int length,double (*CDF)(double))`

Calculate the P-value of Kolmogorov-Smirnov Test of the random sequence and the expected distribution, `CDF` is the function pointer of Cumulative Distribution Function of the expected distribution.
计算随机序列和期望分布的Kolmogorov-Smirnov检验的P-value, `CDF`是期望分布的分布函数的函数指针.


