#pragma once

#include <math_constants.h>
#include <cuda_runtime.h>

//
// Math prototypes
//

template<class REAL> __inline__ __host__ __device__ void gad_sincos( REAL angle, REAL *a, REAL *b );
template<class REAL> __inline__ __host__ __device__ REAL gad_rsqrt( REAL val );


//
// Implementation
//

template<> __inline__ __host__ __device__ void gad_sincos<float>( float angle, float *a, float *b ){ sincosf(angle, a,b); }
template<> __inline__ __host__ __device__ void gad_sincos<double>( double angle, double *a, double *b ){ sincos(angle, a,b); }

template<> __inline__ __host__ __device__ float gad_rsqrt<float>( float val ){ return rsqrtf(val); }
template<> __inline__ __host__ __device__ double gad_rsqrt<double>( double val ){ return rsqrt(val); }
