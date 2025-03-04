/**
 * Description:
 * 
 * Head only pseudorandom engine base on PCG-XSH-RR, 
 * high performance, statistically good and easy to use.
 * 
 * Functions:
 * 
 * Set the seed of PCG random engine and initlize 
 * void PCG32SetSeed(PCG32Struct* status,long long unsigned int seed);
 * 
 * Generate a unsigned type random number in range of [0,0xFFFFFFFF(4294967295)]
 * unsigned PCG32(PCG32Struct* status);
 * 
 * Generate a unsigned type random number that obey uniform distrubution in range of [min,max]
 * unsigned PCG32Uniform(PCG32Struct* status,unsigned min,unsigned max);
 * 
 * Generate a double type random number that obey uniform distrubution in range of [min,max]
 * double PCG32UniformReal(PCG32Struct* status,double min,double max);
 * 
 * Generate a double type random number that obey standard normal distrubution
 * double PCG32StandardNormal(PCG32Struct* status);
 * 
 * Example:
 * 
 * PCG32Struct PCGStatus;
 * PCG32SetSeed(&PCGStatus,time(NULL));
 * double random=PCG32UniformReal(&PCGStatus,-1,999);
 * 
 * Note:
 * 
 * Remember to set the seed.
 * Use its own PCG32Struct in each thread function and set different seed. 
 * 
*/

#ifndef __PCG32_H__
#define __PCG32_H__

#ifdef __cplusplus
extern "C"{
#endif

#include <math.h>

static long long unsigned int const multiplier=6364136223846793005LLU;
static long long unsigned int const increment=1442695040888963407LLU;
static long long unsigned int const PCG32Max=0x00000000FFFFFFFFLLU;
static const unsigned char PCG32False=0;
static const unsigned char PCG32True=1;
static const double PCG32RealScale=(double)1/(double)PCG32Max;

typedef struct _PCG32Struct{
	long long unsigned int state;
	long long unsigned int seed;
	double normalDistributionSaved;
	unsigned normalDistributionSavedValid;
}PCG32Struct;

static inline unsigned rotr32(unsigned x,unsigned r){
	return x>>r|x<<(-r&31);
}

static inline unsigned PCG32(PCG32Struct* status){
	long long unsigned int x=status->state;
	unsigned count=(unsigned)(x>>59);
	status->state=x*multiplier+increment;
	x=x^(x>>18);
	return rotr32((unsigned)(x>>27),count);
}

static inline unsigned PCG32Uniform(PCG32Struct* status,unsigned min,unsigned max){
	if(min>max){
		unsigned tempory=max;
		max=min;
		min=tempory;
	}
	long long unsigned int gap=max-min+1;
	if((gap&(gap-1))==0){
		return (PCG32(status)&(gap-1))+min;
	}
	unsigned range=(unsigned)(((PCG32Max+1)/gap)*gap);
	unsigned random=PCG32(status);
	while(random>range){
		random=PCG32(status);
	}
	return (random%gap)+min;
}

static inline double PCG32UniformReal(PCG32Struct* status,double min,double max){
	return min+((double)PCG32(status))*PCG32RealScale*(max-min);
}

static inline double PCG32StandardNormal(PCG32Struct* status){
	if(status->normalDistributionSavedValid){
		status->normalDistributionSavedValid=PCG32False;
		return status->normalDistributionSaved;
	}
	double u1,u2,S;
	do{
		u1=(double)(PCG32(status))/(double)PCG32Max*2.0-1.0;
		u2=(double)(PCG32(status))/(double)PCG32Max*2.0-1.0;
		S=u1*u1+u2*u2;
	}while(S>1.0||S==0.0);
	const double toMultiple=sqrt(-2.0*log(S)/S);
	status->normalDistributionSavedValid=PCG32True;
	status->normalDistributionSaved=toMultiple*u2;
	return toMultiple*u1;
}

static inline void PCG32Init(PCG32Struct* status){
	status->state=status->seed+increment;
	status->normalDistributionSavedValid=PCG32False;
	PCG32(status);
}

static inline void PCG32SetSeed(PCG32Struct* status,long long unsigned int seed){
	status->seed=seed;
	PCG32Init(status);
}

#ifdef __cplusplus
}
#endif

#endif