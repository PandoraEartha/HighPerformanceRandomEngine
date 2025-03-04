#ifndef __PCG32TEST_H__
#define __PCG32TEST_H__

#ifdef __cplusplus
extern "C"{
#endif

#include "PCG32.h"
#include <stdlib.h>

static const double sqrt2reciprocal=0.70710678118654752440084436210484903928483593768847;

static inline double normalCDF(double x){
	return 0.5*(1+erf(x*sqrt2reciprocal));
}

static inline double uniformCDF(double x){
	return x;
}

static inline int doubleCompare(const void* a,const void* b){
	double aDouble=*((const double*)a);
	double bDouble=*((const double*)b);
	return (aDouble>bDouble)-(aDouble<bDouble);
}

static inline double KSTestPValue(double D,long long unsigned int length){
	double x2=D*D*(double)length;
	return 2*(exp(-2*x2)-exp(-2*4*x2)+exp(-2*9*x2)-exp(-2*16*x2));
}

#define KolmogorovSmirnovTest KSTest

static inline double KSTest(double* data,long long unsigned int length,double (*CDF)(double)){
	qsort(data,length,sizeof(double),doubleCompare);
	double D=0.0;
	for(long long unsigned int index=0;index<length;index=index+1){
		double observationCDF=(double)(index+1)/(double)length;
		double cdf=CDF(data[index]);
		double distance=fabs(observationCDF-cdf);
		if(D<distance){
			D=distance;
		}
	}
	return KSTestPValue(D,length);
}

#ifdef __cplusplus
}
#endif

#endif