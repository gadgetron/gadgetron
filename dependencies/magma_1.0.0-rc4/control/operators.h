/**
 *
 *  @file operators.h
 *
 *  MAGMA (version 1.0) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  November 2010
 *
 **/
#ifndef _MAGMA_OPERATORS_H_
#define _MAGMA_OPERATORS_H_

/*************************************************************
 *              cuDoubleComplex
 */

__host__ __device__ static __inline__ cuDoubleComplex 
operator-(const cuDoubleComplex &a)
{
    return make_cuDoubleComplex(-a.x, -a.y);
}

__host__ __device__ static __inline__ cuDoubleComplex 
operator+(const cuDoubleComplex a, const cuDoubleComplex b)
{
    return make_cuDoubleComplex(a.x + b.x, a.y + b.y);
}

__host__ __device__ static __inline__ void
operator+=(cuDoubleComplex &a, const cuDoubleComplex b)
{
    a.x += b.x; a.y += b.y;
}

__host__ __device__ static __inline__ cuDoubleComplex 
operator-(const cuDoubleComplex a, const cuDoubleComplex b)
{
    return make_cuDoubleComplex(a.x - b.x, a.y - b.y);
}

__host__ __device__ static __inline__ void
operator-=(cuDoubleComplex &a, const cuDoubleComplex b)
{
    a.x -= b.x; a.y -= b.y;
}

__host__ __device__ static __inline__ cuDoubleComplex 
operator*(const cuDoubleComplex a, const cuDoubleComplex b)
{
    return make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y);
}

__host__ __device__ static __inline__ cuDoubleComplex 
operator*(const cuDoubleComplex a, const double s)
{
    return make_cuDoubleComplex(a.x * s, a.y * s);
}

__host__ __device__ static __inline__ cuDoubleComplex 
operator*(const double s, const cuDoubleComplex a)
{
    return make_cuDoubleComplex(a.x * s, a.y * s);
}

__host__ __device__ static __inline__ void 
operator*=(cuDoubleComplex &a, const cuDoubleComplex b)
{
  double tmp = a.y * b.x + a.x * b.y;
  a.x = a.x * b.x - a.y * b.y;
  a.y = tmp;
}

__host__ __device__ static __inline__ void 
operator*=(cuDoubleComplex &a, const double s)
{
    a.x *= s; a.y *= s;
}

/*************************************************************
 *              cuFloatComplex
 */

__host__ __device__ static __inline__ cuFloatComplex 
operator-(const cuFloatComplex &a)
{
    return make_cuFloatComplex(-a.x, -a.y);
}

__host__ __device__ static __inline__ cuFloatComplex 
operator+(const cuFloatComplex a, const cuFloatComplex b)
{
    return make_cuFloatComplex(a.x + b.x, a.y + b.y);
}

__host__ __device__ static __inline__ void
operator+=(cuFloatComplex &a, const cuFloatComplex b)
{
    a.x += b.x; a.y += b.y;
}

__host__ __device__ static __inline__ cuFloatComplex 
operator-(const cuFloatComplex a, const cuFloatComplex b)
{
    return make_cuFloatComplex(a.x - b.x, a.y - b.y);
}

__host__ __device__ static __inline__ void
operator-=(cuFloatComplex &a, const cuFloatComplex b)
{
    a.x -= b.x; a.y -= b.y;
}

__host__ __device__ static __inline__ cuFloatComplex 
operator*(const cuFloatComplex a, const cuFloatComplex b)
{
    return make_cuFloatComplex(a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y);
}

__host__ __device__ static __inline__ cuFloatComplex 
operator*(const cuFloatComplex a, const float s)
{
    return make_cuFloatComplex(a.x * s, a.y * s);
}

__host__ __device__ static __inline__ cuFloatComplex 
operator*(const float s, const cuFloatComplex a)
{
    return make_cuFloatComplex(a.x * s, a.y * s);
}

__host__ __device__ static __inline__ void 
operator*=(cuFloatComplex &a, const cuFloatComplex b)
{
  float tmp = a.y * b.x + a.x * b.y;
  a.x = a.x * b.x - a.y * b.y;
  a.y = tmp;
}

__host__ __device__ static __inline__ void 
operator*=(cuFloatComplex &a, const float s)
{
    a.x *= s; a.y *= s;
}

#endif



