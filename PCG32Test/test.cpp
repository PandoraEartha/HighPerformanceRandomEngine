#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "PCG32.h"
#include <sys/time.h>
#include <string.h>
#include <random>
#include "PCG32Test.h"

const unsigned ShortMax=0xFFFF;
const long long unsigned int TimeTestCase=123456789LLU*33;

const bool SpeedTestOnly=true;

void testUniformReal(PCG32Struct* PCGStatus){
	if(!SpeedTestOnly){
		const unsigned rangeTestCase=1<<30;
		for(long long unsigned int index=0;index<rangeTestCase;index=index+1){
			double min=PCG32UniformReal(PCGStatus,-1,1)*(PCG32(PCGStatus)&ShortMax);
			double max=PCG32UniformReal(PCGStatus,-1,1)*(PCG32(PCGStatus)&ShortMax);
			double random=PCG32UniformReal(PCGStatus,min,max);
			if(min>max){
				double tempory=max;
				max=min;
				min=tempory;
			}
			if(random<min||random>max){
				printf("%lf not in [%lf,%lf]\n",random,min,max);
				return;
			}
		}
	}
	double* uniform=(double*)malloc(sizeof(double)*TimeTestCase);
	struct timeval start,end;
	gettimeofday(&start,NULL);
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		uniform[index]=PCG32UniformReal(PCGStatus,0,1);
	}
	gettimeofday(&end,NULL);
	printf("PCG32UniformReal takes %lf millisecond\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
	if(!SpeedTestOnly){
		printf("PCG32UniformReal KS-Test P-value: %lf\n",KolmogorovSmirnovTest(uniform,TimeTestCase,uniformCDF));
	}

	std::default_random_engine generator(time(NULL));
	std::uniform_real_distribution<double> realDistribution(0,1);
	gettimeofday(&start,NULL);
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		uniform[index]=realDistribution(generator);
	}
	gettimeofday(&end,NULL);
	printf("std::uniform_real_distribution takes %lf millisecond\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
	if(!SpeedTestOnly){
		printf("std::uniform_real_distribution KS-Test P-value: %lf\n",KolmogorovSmirnovTest(uniform,TimeTestCase,uniformCDF));
	}

	free(uniform);
	printf("\n");
}

void testStandardNormal(PCG32Struct* PCGStatus){
	double* normal=(double*)malloc(sizeof(double)*TimeTestCase);
	const double mu=0;
	const double sigma=1;
	struct timeval start,end;
	gettimeofday(&start,NULL);
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		normal[index]=PCG32StandardNormal(PCGStatus);
		normal[index]=normal[index]*sigma+mu;
	}
	gettimeofday(&end,NULL);
	printf("PCG32StandardNormal takes %lf milliseconds\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
	if(!SpeedTestOnly){
		printf("PCG32StandardNormal KS-Test P-value: %lf\n",KSTest(normal,TimeTestCase,normalCDF));
	}

	std::default_random_engine generator(time(NULL));
	std::normal_distribution<double> normalDistribution(mu,sigma);
	gettimeofday(&start,NULL);
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		normal[index]=normalDistribution(generator);
	}
	gettimeofday(&end,NULL);
	printf("std::normal_distribution takes %lf milliseconds\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
	if(!SpeedTestOnly){
		printf("std::normal_distribution KS-Test P-value: %lf\n",KSTest(normal,TimeTestCase,normalCDF));
	}
	free(normal);
	printf("\n");
}

void testUniform(PCG32Struct* PCGStatus){
	if(!SpeedTestOnly){
		const unsigned rangeTestCase=1<<30;
		for(long long unsigned int index=0;index<rangeTestCase;index=index+1){
			unsigned min=PCG32Uniform(PCGStatus,0,77755);
			unsigned max=PCG32Uniform(PCGStatus,999,100033);
			unsigned random=PCG32Uniform(PCGStatus,min,max);
			if(min>max){
				unsigned tempory=max;
				max=min;
				min=tempory;
			}
			if(random<min||random>max){
				printf("%u not in [%u,%u]\n",random,min,max);
				return;
			}
		}
	}
	struct timeval start,end;
	gettimeofday(&start,NULL);
	unsigned PCG32Rand;
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		PCG32Rand=PCG32Uniform(PCGStatus,0,10000-1);
	}
	gettimeofday(&end,NULL);
	printf("PCG32Uniform takes %lf milliseconds\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);

	std::default_random_engine generator(time(NULL));
	gettimeofday(&start,NULL);
	unsigned stdRand;
	std::uniform_int_distribution<unsigned> unifromDistribution(0,10000-1);
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		stdRand=unifromDistribution(generator);
	}
	gettimeofday(&end,NULL);
	printf("std::uniform_int_distribution takes %lf milliseconds\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
	printf("\n");
}

void testBase(PCG32Struct* PCGStatus){
	struct timeval start,end;
	gettimeofday(&start,NULL);
	unsigned PCG32Rand;
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		PCG32Rand=PCG32(PCGStatus);
	}
	gettimeofday(&end,NULL);
	printf("PCG32 takes %lf milliseconds\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);

	std::default_random_engine generator(time(NULL));
	gettimeofday(&start,NULL);
	unsigned stdRand;
	for(long long unsigned int index=0;index<TimeTestCase;index=index+1){
		stdRand=generator();
	}
	gettimeofday(&end,NULL);
	printf("std::default_random_engine takes %lf milliseconds\n",(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
	printf("\n");
}

int main(int argc, char const *argv[]){
	PCG32Struct PCGStatus;
	PCG32SetSeed(&PCGStatus,time(NULL));
	for(unsigned index=0;index<10;index=index+1){
		testUniformReal(&PCGStatus);
		PCG32SetSeed(&PCGStatus,time(NULL));
		testStandardNormal(&PCGStatus);
		PCG32SetSeed(&PCGStatus,time(NULL));
		testUniform(&PCGStatus);
		PCG32SetSeed(&PCGStatus,time(NULL));
		testBase(&PCGStatus);
	}
	return 0;
}
