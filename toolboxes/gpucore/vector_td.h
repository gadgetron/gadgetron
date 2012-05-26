#pragma once

#include "gpucore_defines.h"

template< class T, unsigned int D > class vector_td
{
public:

  T vec[D];

  __inline__ __gad_host__ __gad_device__ T& operator[](const int i){
    return vec[i];
  }

  __inline__ __gad_host__ __gad_device__ const T& operator[](const int i) const {
    return vec[i];
  }
};

//
// Some typedefs for convenience (templated typedefs are not (yet) available in C++)
//

template< class REAL, unsigned int D > struct reald{
  typedef vector_td< REAL, D > Type;
};

template< unsigned int D > struct intd{
  typedef vector_td< int, D > Type;
};

template< unsigned int D > struct uintd{
  typedef vector_td< unsigned int, D > Type;
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

  __inline__ __gad_host__ __gad_device__ vector_td(){}

  __inline__ __gad_host__ __gad_device__ vector_td(T x){
    vec[0]=x;
  }

  __inline__ __gad_host__ __gad_device__ T& operator[](const int i){
    return vec[i];
  }

  __inline__ __gad_host__ __gad_device__ const T& operator[](const int i) const {
    return vec[i];
  }
};

template<class T> class vector_td<T,2>
{
public:

  T vec[2];

  __inline__ __gad_host__ __gad_device__ vector_td(){}

  __inline__ __gad_host__ __gad_device__ vector_td(T x, T y){
    vec[0]=x;
    vec[1]=y;
  }

  __inline__ __gad_host__ __gad_device__ T& operator[](const int i){
    return vec[i];
  }

  __inline__ __gad_host__ __gad_device__ const T& operator[](const int i) const {
    return vec[i];
  }
};

template<class T> class vector_td<T,3>
{
public:

  T vec[3];

  __inline__ __gad_host__ __gad_device__ vector_td(){}

  __inline__ __gad_host__ __gad_device__ vector_td(T x, T y,T z){
    vec[0]=x;
    vec[1]=y;
    vec[2]=z;
  }

  __inline__ __gad_host__ __gad_device__ T& operator[](const int i){
    return vec[i];
  }

  __inline__ __gad_host__ __gad_device__ const T& operator[](const int i) const {
    return vec[i];
  }
};

template<class T> class vector_td<T,4>
{
public:

  T vec[4];

  __inline__ __gad_host__ __gad_device__ vector_td(){}

  __inline__ __gad_host__ __gad_device__ vector_td(T x, T y,T z,T w){
    vec[0]=x;
    vec[1]=y;
    vec[2]=z;
    vec[3]=w;
  }

  __inline__ __gad_host__ __gad_device__ T& operator[](const int i){
    return vec[i];
  }

  __inline__ __gad_host__ __gad_device__ const T& operator[](const int i) const {
    return vec[i];
  }
};

typedef vector_td<double,1> doubled1;
typedef vector_td<double,2> doubled2;
typedef vector_td<double,3> doubled3;
typedef vector_td<double,4> doubled4;

typedef vector_td<float,1> floatd1;
typedef vector_td<float,2> floatd2;
typedef vector_td<float,3> floatd3;
typedef vector_td<float,4> floatd4;

typedef vector_td<int,1> intd1;
typedef vector_td<int,2> intd2;
typedef vector_td<int,3> intd3;
typedef vector_td<int,4> intd4;

typedef vector_td<unsigned int,1> uintd1;
typedef vector_td<unsigned int,2> uintd2;
typedef vector_td<unsigned int,3> uintd3;
typedef vector_td<unsigned int,4> uintd4;
