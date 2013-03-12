#pragma once

#if defined(GTCUDA_FOUND) || defined(__host__) || defined(__CUDACC__)
#include <host_defines.h>
#else
//#define __inline__ inline
//#define __no_return__
//#define __noinline__
//#define __forceinline__
//#define __align__(n)
#define __thread__
#define __host__
#define __device__
#define __global__
#define __shared__
#define __constant__

#endif
