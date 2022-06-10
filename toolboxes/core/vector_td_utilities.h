/** \file vector_td_utilities.h
    \brief The class vector_td defines a D-dimensional vector of type T.

    The class vector_td defines a D-dimensional vector of type T.
    It is used in the Gadgetron to represent small (one- to four-dimensional) vectors only.
    For larger vectors consider using the NDArray class instead.
    The vector_td class can be used on both the cpu and gpu.
    The accompanying headers vector_td_opeators.h and vector_td_utilities.h define most of the functionality.
*/

#pragma once

#include "vector_td.h"
#include "vector_td_operators.h"
#include "real_utilities.h"
#include "core_defines.h"
#include "complext.h"

#include <float.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#ifndef __CUDA_ARCH__ // workaround for nvcc
using std::ceil;  
using std::floor; 
using std::abs;   
using std::sqrt;
#endif

namespace Gadgetron{

  // Windows/Cuda has some issues when using min and max.
  // For now we define our own implementation

  template <class T> __inline__ __host__ __device__ const T& _vector_td_min (const T& a, const T& b) {
    return (a>b)?b:a;
  }
  template <class T> __inline__ __host__ __device__ const T& _vector_td_max (const T& a, const T& b) {
    return (a<b)?b:a;
  }

  //
  // In-place operations
  //

  template<class T, unsigned int D> __inline__ __host__ __device__
  void clear( vector_td<T,D> &vec, const T &val = T(0) )
  {
    for (unsigned int i=0; i<D; i++) {
      vec[i] = val;
    }
  }
  
  //
  // Component-wise math operations
  //

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<T,D> abs( const vector_td<T,D>& vec )
  {
    vector_td<T,D> res;
    for (unsigned int i=0; i<D; i++) {
      res[i] = ::abs(vec[i]);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<int,D> sgn( const vector_td<T,D>& vec )
  {
    vector_td<int,D> res;
    for (unsigned int i=0; i<D; i++) {
      res[i] = sgn(vec[i]);
    }
    return res;
  }

  template<class REAL, unsigned int D> __inline__ __host__ __device__
  vector_td<REAL,D> ceil( const vector_td<REAL,D>& vec )
  {
    vector_td<REAL,D> res;
    for (unsigned int i=0; i<D; i++) {
      res[i] = ::ceil(vec[i]);
    }
    return res;
  }

  template<class REAL, unsigned int D> __inline__ __host__ __device__
  vector_td<REAL,D> floor( const vector_td<REAL,D>& vec )
  {
    vector_td<REAL,D> res;
    for (unsigned int i=0; i<D; i++) {
      res[i] = ::floor(vec[i]);
    }
    return res;
  }


  //
  // Grid <-> index transformations
  //

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<T,D> idx_to_co( T idx, const vector_td<T,D> dims )
  {
    vector_td<T,D> co;
    auto idx_tmp = idx;
    for (int i=0; i<D; i++) {
      co[i] = idx_tmp%dims[i];
      idx_tmp -= co[i];
      idx_tmp /= dims[i];
    }
    return co;
  } 


  template<class T, unsigned int D> __inline__ __host__ __device__
  T co_to_idx( const vector_td<T,D> co, const vector_td<T,D> dims )
  {
    T idx = 0;
    T block_size = 1;
    for (int i=0; i<D; i++) {
      idx += (block_size*co[i]);
      block_size *= dims[i];
    }
    return idx;
  }

  
  template<unsigned int D> __inline__ __host__ __device__
  unsigned int co_to_idx( const vector_td<unsigned int,D> co, 
                          const vector_td<unsigned int,D> dims, 
                          const vector_td<unsigned int,D> order )
  {
    unsigned int idx = 0;
    unsigned int block_size = 1;
    for (unsigned int i=0; i<D; i++){
      idx += (block_size*co.d[order[i]]);
      block_size *= dims.d[order[i]];
    }
    return idx;
  } 

  template<unsigned int D> __inline__ __host__ __device__
  size_t co_to_idx( const vector_td<size_t,D> co, 
                    const vector_td<size_t,D> dims, 
                    const vector_td<unsigned int,D> order )
  {
    size_t idx = 0;
    size_t block_size = 1;
    for (unsigned int i=0; i<D; i++){
      idx += (block_size*co.d[order[i]]);
      block_size *= dims.d[order[i]];
    }
    return idx;
  } 

  template<int D> __inline__ __host__ __device__ 
  int co_to_idx( const vector_td<int,D> co, 
                 const vector_td<int,D> dims, 
                 const vector_td<unsigned int,D> order )
  {
    int idx = 0;
    int block_size = 1;
    for (unsigned int i=0; i<D; i++){
      idx += (block_size*co.d[order[i]]);
      block_size *= dims.d[order[i]];
    }
    return idx;
  } 

  template<unsigned int D> __inline__ __host__ __device__
  long long co_to_idx( const vector_td<long long,D> co, 
                       const vector_td<long long,D> dims, 
                       const vector_td<unsigned int,D> order )
  {
    long long idx = 0;
    long long block_size = 1;
    for (unsigned int i=0; i<D; i++){
      idx += (block_size*co.d[order[i]]);
      block_size *= dims.d[order[i]];
    }
    return idx;
  } 

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<T,D> counting_vec()
  {
    vector_td<T,D> res;
    for(unsigned int i=0; i<D; i++) {
      res[i]=T(i);
    }
    return res;
  }

  //
  // Conversion between vector_td and std::vector
  //

  template<class T, unsigned int D> inline
  std::vector<T> to_std_vector(const vector_td<T,D>& vec )
  {
    std::vector<T> out(D);
    for(unsigned int i=0; i<D; i++ )
      out[i] = vec[i];
    return out;
  }

  template<class T, unsigned int D> inline
  vector_td<T,D> from_std_vector(const std::vector<T> &_vector )
  {
    vector_td<T,D> out;
    for( unsigned int i=0; i<D; i++ ){
      if( i<_vector.size() )
        out[i] = _vector[i];
      else
        out[i] = T(1);
    }
    return out;
  }

  //
  // Reductions on vector_td<T,D>
  //

  template<class T, unsigned int D> __inline__ __host__ __device__
  T prod( const vector_td<T,D>& vec )
  {
    T res = vec[0];
    for (unsigned int i=1; i<D; i++){
      res *= vec[i];
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  T sum( const vector_td<T,D>& vec )
  {
    T res = vec[0];
    for (unsigned int i=1; i<D; i++){
      res += vec[i];
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  T dot( const vector_td<T,D>& vec1, const vector_td<T,D>& vec2 )
  {
    T res = (vec1[0]*vec2[0]);
    for (unsigned int i=1; i<D; i++){
      res += (vec1[i]*vec2[i]);
    }
    return res;
  }

  template<class REAL, unsigned int D> __inline__ __host__ __device__
  complext<REAL> dot(const vector_td<complext<REAL>, D>& vec1, const vector_td<REAL, D>& vec2)
  {
	  complext<REAL> res = (vec1[0] * vec2[0]);
	  for (unsigned int i = 1; i<D; i++){
		  res += (vec1[i] * vec2[i]);
	  }
	  return res;

  }

  template<class REAL, unsigned int D> __inline__ __host__ __device__
  complext<REAL> dot(const vector_td<REAL, D>& vec1, const vector_td<complext<REAL>, D>& vec2)
  {
	  complext<REAL> res = (vec1[0] * vec2[0]);
	  for (unsigned int i = 1; i<D; i++){
		  res += (vec1[i] * vec2[i]);
	  }
	  return res;

  }
  template<class T, unsigned int D> __inline__ __host__ __device__
  T max( const vector_td<T,D>& vec )
  {
    T res = vec[0];
    for (unsigned int i=1; i<D; i++){
      res = _vector_td_max(res,vec[i]);
    }
    return res;
  }
  
  template<class T, unsigned int D> __inline__ __host__ __device__
  T min( const vector_td<T,D>& vec )
  {
    T res = vec[0];
    for (unsigned int i=1; i<D; i++){
      res = _vector_td_min(res,vec[i]);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<T,D> amin( const vector_td<T,D>& vec1, const vector_td<T,D>& vec2)
  {
    vector_td<T,D> res;
    for (unsigned int i=0; i<D; i++){
      res[i] = _vector_td_min(vec1[i],vec2[i]);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<T,D> amax( const vector_td<T,D>& vec1, const vector_td<T,D>& vec2)
  {
    vector_td<T,D> res;
    for (unsigned int i=0; i<D; i++){
      res[i] = _vector_td_max(vec1[i],vec2[i]);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<T,D> amin( const vector_td<T,D>& vec1, T val)
  {
    vector_td<T,D> res;
    for (unsigned int i=0; i<D; i++){
      res[i] = _vector_td_min(vec1[i],val);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  vector_td<T,D> amax( const vector_td<T,D>& vec1, T val )
  {
    vector_td<T,D> res;
    for (unsigned int i=0; i<D; i++){
      res[i] = _vector_td_max(vec1[i],val);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  T max_not_nan( const vector_td<T,D>& vec )
  {
    unsigned int i=0;
    while (isnan(vec[i])) i++;
    if (i >= D) return 0;
    T res = vec[i];
    for (++i; i<D; i++){
      if (!isnan(vec[i])) res = _vector_td_max(res,vec[i]);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  T min_not_nan( const vector_td<T,D>& vec )
  {
    unsigned int i=0;
    while (isnan(vec[i])) i++;
    T res = vec[i];
    for (++i; i<D; i++){
      if (!isnan(vec[i])) res = _vector_td_min(res,vec[i]);
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  unsigned int argmin( const vector_td<T,D>& vec )
  {
    unsigned int res= 0;
    for (unsigned int i=1; i<D; i++){
      if (vec[i] < vec[res] ) res = i;
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  unsigned int argmin_not_nan( const vector_td<T,D>& vec )
  {
    unsigned int res= 0;
    for (unsigned int i=1; i<D; i++){
      if (vec[i] < vec[res] && !isnan(vec[i])) res = i;
    }
    return res;
  }

  template<class T, unsigned int D> __inline__ __host__ __device__
  unsigned int argmax( const vector_td<T,D>& vec )
  {
    unsigned int res= 0;
    for (unsigned int i=1; i<D; i++){
      if (vec[i] > vec[res] ) res = i;
    }
    return res;
  }

  //
  // Reductions on reald<REAL,D>
  //

  template<class REAL, unsigned int D> __inline__ __host__ __device__
  typename realType<REAL>::Type norm_squared( const vector_td<REAL,D> vec )
  {
    typename realType<REAL>::Type res(0);
    for (unsigned int i=0; i<D; i++){
      res += norm(vec[i]);
    }
    return res;
  }

  template<class REAL, unsigned int D> __inline__ __host__ __device__
  typename realType<REAL>::Type norm( const vector_td<REAL,D> vec )
  {
    return ::sqrt(norm_squared<REAL,D>(vec));
  }

}
