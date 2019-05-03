/** \file real_utilities.h
    \brief A simple template based interface to some common C float/double constants to ease writing of templated code.
*/

#pragma once

#include "core_defines.h"

#ifdef _USE_MATH_DEFINES
#include <math.h>
#else
#define _USE_MATH_DEFINES
#include <math.h>
#undef _USE_MATH_DEFINES
#endif

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

template<> __inline__ __host__ __device__ float get_pi(){ return (float)M_PI; }
template<> __inline__ __host__ __device__ double get_pi(){ return M_PI; }

