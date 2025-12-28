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
 * Generate a double type random number that obey uniform distrubution in range of [min,max)
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

#if defined(__CUDACC__)||defined(__CUDA_ARCH__)||defined(__CUDA_LIBDEVICE__)
    #define PCG32_CUDA 1
#else
    #define PCG32_CUDA 0
#endif

#if PCG32_CUDA
    #define PCG32_HOST_DEVICE __host__ __device__
    #define PCG32_DEVICE               __device__
    #define PCG32_HOST        __host__
	#include <cuda_runtime_api.h>
#else
    #define PCG32_HOST_DEVICE
    #define PCG32_DEVICE
    #define PCG32_HOST
	#include <math.h>
	#include <stdbool.h>
#endif

#if PCG32_CUDA
	#define multiplier 6364136223846793005LLU
	#define increment  1442695040888963407LLU
	#define PCG32Max   0x00000000FFFFFFFFLLU
	#define PCG32False 0
	#define PCG32True  1
	#define PCG32RealScale 1.0/(PCG32Max+1)
#else
	static long long unsigned int const multiplier=6364136223846793005LLU;
	static long long unsigned int const increment=1442695040888963407LLU;
	static long long unsigned int const PCG32Max=0x00000000FFFFFFFFLLU;
	static const unsigned char PCG32False=0;
	static const unsigned char PCG32True=1;
	static const double PCG32RealScale=(double)1/(double)(PCG32Max+1);
#endif

#ifdef __cplusplus
extern "C"{
#endif



typedef struct _PCG32Struct{
	long long unsigned int state;
	long long unsigned int seed;
	double normalDistributionSaved;
	double gammaD;
	double gammaC;
	double gammaBeta;
	double gammaNegativeSqrt9AlphaSub3;
	unsigned normalDistributionSavedValid;
}PCG32Struct;

PCG32_HOST_DEVICE static inline unsigned rotr32(unsigned x,unsigned r){
	return x>>r|x<<(-r&31);
}

PCG32_HOST_DEVICE static inline unsigned PCG32(PCG32Struct* status){
	long long unsigned int x=status->state;
	unsigned count=(unsigned)(x>>59);
	status->state=x*multiplier+increment;
	x=x^(x>>18);
	return rotr32((unsigned)(x>>27),count);
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
	unsigned range=(unsigned)(((PCG32Max+1)/gap)*gap);
	unsigned random=PCG32(status);
	while(random>range){
		random=PCG32(status);
	}
	return (random%gap)+min;
}

// max can not smaller than min and gap can not be power of 2
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_Strict(PCG32Struct* status,const unsigned min,const unsigned max){
	unsigned gap=max-min+1;
	unsigned range=(unsigned)(((PCG32Max+1)/gap)*gap);
	unsigned random=PCG32(status);
	while(random>range){
		random=PCG32(status);
	}
	return (random%gap)+min;
}

// max can not smaller than min
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_MaxBiggerThanMin(PCG32Struct* status,const unsigned min,const unsigned max){
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

PCG32_HOST_DEVICE static inline double PCG32UniformReal(PCG32Struct* status,double min,double max){
	return min+((double)PCG32(status))*PCG32RealScale*(max-min);
}

PCG32_HOST_DEVICE static inline double PCG32StandardNormal(PCG32Struct* status){
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

// Currently, only the algorithm for a >= 1 has been implemented
PCG32_HOST_DEVICE static inline bool PCG32GammaInit(PCG32Struct* status,double alpha,double beta){
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

PCG32_HOST_DEVICE static inline void PCG32Init(PCG32Struct* status){
	status->state=status->seed+increment;
	status->normalDistributionSavedValid=PCG32False;
	PCG32(status);
}

PCG32_HOST_DEVICE static inline void PCG32SetSeed(PCG32Struct* status,long long unsigned int seed){
	status->seed=seed;
	PCG32Init(status);
}

#ifdef __cplusplus
}
#endif

#if defined(__cplusplus)||PCG32_CUDA

template<typename Type>
PCG32_HOST_DEVICE static inline void PCG32SWAP(Type* array,long long unsigned int index0,long long unsigned int index1){
	const Type tempory=array[index0];                                                                            
    array[index0]=array[index1];                                                                                 
    array[index1]=tempory;  
}

template<typename Type>
PCG32_HOST_DEVICE static inline void PCG32UniformShuffle(PCG32Struct* status,Type* array,long long unsigned int length){
	if(length>1){                                                                                
    	for(long long unsigned int index=0;index<length-1;index=index+1){                        
            PCG32SWAP(array,index,PCG32Uniform_MaxBiggerThanMin(status,index,length-1));
        }                                                                                        
    }   
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
        long long*:			 PCG32SWAP_long_long,         \
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

#endif

#endif