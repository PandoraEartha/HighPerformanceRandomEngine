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

#define PCG32_INT_MAX      0x7FFFFFFFU
#define PCG32_UNSIGNED_MAX 0xFFFFFFFFU

#define PCG32BINOMIAL_SMALLMEAN    14
#define PCG32BINOMIAL_MAXITERATION 110
#define PCG32BINOMIAL_FARFROMMENA  20

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
	unsigned uniformStrictGap;
	unsigned uniformStrictRange;
	unsigned uniformStrictMin;
}PCG32Struct;

PCG32_HOST_DEVICE static inline unsigned rotr32(unsigned x,unsigned r){
	return x>>r|x<<(-r&31);
}

PCG32_HOST_DEVICE static inline unsigned PCG32(PCG32Struct* status){
	long long unsigned int x=status->state;
	unsigned count=(unsigned)(x>>59);
	status->state=x*PCG32MULTIPLIER+PCG32INCREMENT;
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
	unsigned range=(unsigned)(((PCG32MAX+1)/gap)*gap);
	unsigned random=PCG32(status);
	while(random>range){
		random=PCG32(status);
	}
	return (random%gap)+min;
}

// max can not smaller than min and gap can not be power of 2
PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_Strict(PCG32Struct* status,const unsigned min,const unsigned max){
	unsigned gap=max-min+1;
	unsigned range=(unsigned)(((PCG32MAX+1)/gap)*gap);
	unsigned random=PCG32(status);
	while(random>range){
		random=PCG32(status);
	}
	return (random%gap)+min;
}

PCG32_HOST_DEVICE static inline unsigned PCG32Uniform_StrictRangeUnchanged(PCG32Struct* status,const unsigned min,const unsigned max){
	unsigned random=PCG32(status);
	while(random>status->uniformStrictRange){
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
	while(random>range){
		random=PCG32(status);
	}
	return (random%gap)+min;
}

PCG32_HOST_DEVICE static inline double PCG32UniformReal(PCG32Struct* status,double min,double max){
	return min+((double)PCG32(status))*PCG32REAL_SCALE*(max-min);
}

PCG32_HOST_DEVICE static inline double PCG32StandardNormal(PCG32Struct* status){
	if(status->normalDistributionSavedValid){
		status->normalDistributionSavedValid=PCG32FALSE;
		return status->normalDistributionSaved;
	}
	double u1,u2,S;
	do{
		u1=(double)(PCG32(status))/(double)PCG32MAX*2.0-1.0;
		u2=(double)(PCG32(status))/(double)PCG32MAX*2.0-1.0;
		S=u1*u1+u2*u2;
	}while(S>1.0||S==0.0);
	const double toMultiple=sqrt(-2.0*log(S)/S);
	status->normalDistributionSavedValid=PCG32TRUE;
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
	status->state=status->seed+PCG32INCREMENT;
	status->normalDistributionSavedValid=PCG32FALSE;
	PCG32(status);
}

PCG32_HOST_DEVICE static inline void PCG32SetSeed(PCG32Struct* status,long long unsigned int seed){
	status->seed=seed;
	PCG32Init(status);
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

PCG32_HOST_DEVICE static inline double PCG32Stirling(const double x){
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

PCG32_HOST_DEVICE static inline unsigned PCG32Binomial(PCG32Struct* status,double probability,const unsigned repeatUnsigned){
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
				accept=xm*log(f1/x1)+(repeat-m+0.5)*log(z1/w1)+(success-m)*log(w1*probability/(x1*q))+PCG32Stirling(f1)+PCG32Stirling(z1)-PCG32Stirling(x1)-PCG32Stirling(w1);
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