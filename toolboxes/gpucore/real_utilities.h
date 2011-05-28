#pragma once

#include <math_constants.h>
#include <math_functions.h>
#include <float.h>

//
// Get scalar limits of operation
//

template<class T> __inline__ __host__ __device__ T get_min();
template<class T> __inline__ __host__ __device__ T get_max();
template<class T> __inline__ __host__ __device__ T get_epsilon();

//
// Math prototypes
//

template<class REAL> __inline__ __host__ __device__ void gad_sincos( REAL angle, REAL *a, REAL *b );
template<class REAL> __inline__ __host__ __device__ REAL gad_rsqrt( REAL val );
template<class REAL> __inline__ __device__ REAL get_pi();


//
// Implementation
//

template<> __inline__ __host__ __device__ float get_min<float>()
{
  return FLT_MIN;
}

template<> __inline__ __host__ __device__ double get_min<double>()
{
  return DBL_MIN;
}

template<> __inline__ __host__ __device__ float get_max<float>()
{
  return FLT_MAX;
}

template<> __inline__ __host__ __device__ double get_max<double>()
{
  return DBL_MAX;
}

template<> __inline__ __host__ __device__ float get_epsilon<float>()
{
  return FLT_EPSILON;
}

template<> __inline__ __host__ __device__ double get_epsilon<double>()
{
  return DBL_EPSILON;
}

template<> __inline__ __host__ __device__ void gad_sincos<float>( float angle, float *a, float *b ){ sincosf(angle, a,b); }
template<> __inline__ __host__ __device__ void gad_sincos<double>( double angle, double *a, double *b ){ sincos(angle, a,b); }

template<> __inline__ __host__ __device__ float gad_rsqrt<float>( float val ){ return rsqrtf(val); }
template<> __inline__ __host__ __device__ double gad_rsqrt<double>( double val ){ return rsqrt(val); }

template<> __inline__ __device__ float get_pi(){ return CUDART_PI_F; }
template<> __inline__ __device__ double get_pi(){ return CUDART_PI; }
