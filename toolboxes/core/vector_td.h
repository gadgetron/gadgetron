/** \file vector_td.h
    \brief The class vector_td defines a D-dimensional vector of type T.

    The class vector_td defines a D-dimensional vector of type T.
    It is used in the Gadgetron to represent short vectors.
    I.e. it is purposedly templetated with dimensionality D as type unsigned int instead of size_t.
    For larger vectors consider using the NDArray class instead (or a std::vector).
    The vector_td class can be used on both the cpu and gpu.
    The accompanying headers vector_td_opeators.h and vector_td_utilities.h define most of the functionality.
    Note that vector_td should not be used to represent complex numbers. For that we provide the custom class complext instead.
*/

#pragma once

#include "core_defines.h"

#include <stdlib.h> // for size_t

#ifdef max
#undef max
#endif // max

namespace Gadgetron{

  template<class T, unsigned int D> class vector_td
  {
  public:

    T vec[D];
    __inline__ __host__ __device__ vector_td(){};

#if __cplusplus > 199711L
    template <typename... X>
    constexpr __inline__ __host__ __device__ vector_td(X... xs) : vec{xs...} { }
#endif
    __inline__ __host__ __device__ vector_td(const vector_td & other){
       	for (unsigned int i = 0; i < D; i++)
           	vec[i] = other[i];
        }

    template<class T2> __inline__ __host__ __device__ explicit vector_td(const vector_td<T2,D> & other){
    	for (unsigned int i = 0; i < D; i++)
        	vec[i] = (T) other[i];
     }

    __inline__ __host__ __device__ explicit vector_td(T x){
    	for (unsigned int i = 0; i < D; i++)
               	vec[i] = x;
		 }
    __inline__ __host__ __device__ T& operator[](const unsigned int i)
    {
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const unsigned int i) const
    {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const unsigned int i)
    {
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const unsigned int i) const
    {
      return vec[i];
    }
  };

  //
  // Some typedefs for convenience (templated typedefs are not (yet) available in C++)
  //

  template< class REAL, unsigned int D > struct reald{
    typedef vector_td< REAL, D > Type;
  };

  template< unsigned int D > struct uintd{
    typedef vector_td< unsigned int, D > Type;
  };

  template< unsigned int D > struct uint64d{
    typedef vector_td< size_t, D > Type;
  };

  template< unsigned int D > struct intd{
    typedef vector_td< int, D > Type;
  };

  template< unsigned int D > struct int64d{
    typedef vector_td< long long, D > Type;
  };

  template< unsigned int D > struct floatd{
    typedef typename reald< float, D >::Type Type;
  };

  template< unsigned int D > struct doubled{
    typedef typename reald< double, D >::Type Type;
  };

  template<class T> class vector_td<T,1>
  {
  public:

    T vec[1];

    __inline__ __host__ __device__ vector_td(const vector_td & other){
					vec[0] = other[0];
		 }
    template<class T2> __inline__ __host__ __device__ explicit vector_td(const vector_td<T2,1> & other){
    	vec[0] = (T) other[0];
    }

    __inline__ __host__ __device__ vector_td(){}

    __inline__ __host__ __device__ vector_td(T x){ // Not explicit because we actually want to be able to do implicit conversions here.
      vec[0]=x;
    }

    __inline__ __host__ __device__ T& operator[](const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const unsigned int i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const unsigned int i) const {
      return vec[i];
    }
  };

  template<class T> class vector_td<T,2>
  {
  public:

    T vec[2];

    __inline__ __host__ __device__ vector_td(const vector_td & other){
      	for (unsigned int i = 0; i < 2; i++)
          	vec[i] = other[i];
		 }

    //__inline__ __host__ __device__ explicit vector_td(T (&other)[2]) : vec(other){};

    template<class T2> __inline__ __host__ __device__ explicit vector_td(const vector_td<T2,2> & other){
    	for (unsigned int i = 0; i < 2; i++)
        	vec[i] = (T) other[i];
     }
#if __cplusplus > 199711L

    constexpr __inline__ __host__ __device__ vector_td( T x, T y) : vec{x,y} { }
//    template <typename... X>
//    constexpr __inline__ __host__ __device__ vector_td(X... xs) : vec{xs...} { }
#else
    __inline__ __host__ __device__ vector_td(T x, T y){
    	vec[0] = x;
    	vec[1] = y;
    }
#endif
    __inline__ __host__ __device__ vector_td(){}

    __inline__ __host__ __device__ explicit vector_td(T x){
      vec[0]=x;
      vec[1]=x;
    }
    __inline__ __host__ __device__ T& operator[](const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const unsigned int i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const unsigned int i) const {
      return vec[i];
    }
  };

  template<class T> class vector_td<T,3>
  {
  public:

    T vec[3];

    __inline__ __host__ __device__ vector_td(const vector_td & other){
      	for (unsigned int i = 0; i < 3; i++)
          	vec[i] = other[i];
		 }
    template<class T2> __inline__ __host__ __device__ explicit vector_td(const vector_td<T2,3> & other){
    	for (unsigned int i = 0; i < 3; i++)
        	vec[i] = (T) other[i];
     }
    __inline__ __host__ __device__ vector_td(){}

#if __cplusplus > 199711L
      constexpr __inline__ __host__ __device__ vector_td( T x, T y, T z) : vec{x,y,z} { }
#else
      __inline__ __host__ __device__ vector_td(T x, T y, T z){
    	vec[0] = x;
    	vec[1] = y;
    	vec[2] = z;
    }
#endif

    __inline__ __host__ __device__ explicit vector_td(T x){
      vec[0]=x;
      vec[1]=x;
      vec[2]=x;
    }

    __inline__ __host__ __device__ T& operator[](const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const unsigned int i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const unsigned int i) const {
      return vec[i];
    }
  };

  template<class T> class vector_td<T,4>
  {
  public:

    T vec[4];

    __inline__ __host__ __device__ vector_td(const vector_td & other){
    	for (unsigned int i = 0; i < 4; i++)
        	vec[i] = other[i];
     }
    //__inline__ __host__ __device__ explicit vector_td(T (&other)[4]) : vec(other){};
    template<class T2> __inline__ __host__ __device__ explicit vector_td(const vector_td<T2,4> & other){
    	for (unsigned int i = 0; i < 4; i++)
        	vec[i] = (T) other[i];
     }

#if __cplusplus > 199711L
    constexpr __inline__ __host__ __device__ vector_td( T x, T y, T z, T t) : vec{x,y,z,t} { }
#else
    __inline__ __host__ __device__ vector_td(T x, T y, T z, T t){
    	vec[0] = x;
    	vec[1] = y;
    	vec[2] = z;
    	vec[3] = t;
    }
#endif

    __inline__ __host__ __device__ vector_td(){}

    __inline__ __host__ __device__ explicit vector_td(T x){
      vec[0]=x;
      vec[1]=x;
      vec[2]=x;
      vec[3]=x;
    }

    __inline__ __host__ __device__ T& operator[](const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const unsigned int i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const unsigned int i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const unsigned int i) const {
      return vec[i];
    }
  };

  typedef vector_td<unsigned int,1> uintd1;
  typedef vector_td<unsigned int,2> uintd2;
  typedef vector_td<unsigned int,3> uintd3;
  typedef vector_td<unsigned int,4> uintd4;

  typedef vector_td<size_t,1> uint64d1;
  typedef vector_td<size_t,2> uint64d2;
  typedef vector_td<size_t,3> uint64d3;
  typedef vector_td<size_t,4> uint64d4;

  typedef vector_td<int,1> intd1;
  typedef vector_td<int,2> intd2;
  typedef vector_td<int,3> intd3;
  typedef vector_td<int,4> intd4;

  typedef vector_td<long long,1> int64d1;
  typedef vector_td<long long,2> int64d2;
  typedef vector_td<long long,3> int64d3;
  typedef vector_td<long long,4> int64d4;

  typedef vector_td<float,1> floatd1;
  typedef vector_td<float,2> floatd2;
  typedef vector_td<float,3> floatd3;
  typedef vector_td<float,4> floatd4;

  typedef vector_td<double,1> doubled1;
  typedef vector_td<double,2> doubled2;
  typedef vector_td<double,3> doubled3;
  typedef vector_td<double,4> doubled4;
}
