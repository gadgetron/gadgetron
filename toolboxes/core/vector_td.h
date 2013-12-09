/** \file vector_td.h
    \brief The class vector_td defines a D-dimensional vector of type T.

    The class vector_td defines a D-dimensional vector of type T.
    It is used in the Gadgetron to represent small (one- to four-dimensional) vectors only.
    For larger vectors consider using the NDArray class instead.
    The vector_td class can be used on both the cpu and gpu.
    The accompanying headers vector_td_opeators.h and vector_td_utilities.h define most of the functionality.
*/

#pragma once

#include "core_defines.h"

namespace Gadgetron{

  template< class T, unsigned long long D > class vector_td
  {
  public:

    T vec[D];

    __inline__ __host__ __device__ T& operator[](const long long i)
    {
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const long long i) const
    {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const long long i)
    {
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const long long i) const
    {
      return vec[i];
    }
  };

  //
  // Some typedefs for convenience (templated typedefs are not (yet) available in C++)
  //

  template< class REAL, unsigned long long D > struct reald{
    typedef vector_td< REAL, D > Type;
  };

  template< unsigned long long D > struct intd{
    typedef vector_td< long long, D > Type;
  };

  template< unsigned long long D > struct uintd{
    typedef vector_td< unsigned long long, D > Type;
  };

  template< unsigned long long D > struct floatd{
    typedef typename reald< float, D >::Type Type;
  };

  template< unsigned long long D > struct doubled{
    typedef typename reald< double, D >::Type Type;
  };

  template<class T> class vector_td<T,1>
  {
  public:

    T vec[1];

    __inline__ __host__ __device__ vector_td(){}

    __inline__ __host__ __device__ vector_td(T x){
      vec[0]=x;
    }

    __inline__ __host__ __device__ T& operator[](const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const long long i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const long long i) const {
      return vec[i];
    }
  };

  template<class T> class vector_td<T,2>
  {
  public:

    T vec[2];

    __inline__ __host__ __device__ vector_td(){}

    __inline__ __host__ __device__ vector_td(T x, T y){
      vec[0]=x;
      vec[1]=y;
    }

    __inline__ __host__ __device__ vector_td(T x){
      vec[0]=x;
      vec[1]=x;
    }
    __inline__ __host__ __device__ T& operator[](const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const long long i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const long long i) const {
      return vec[i];
    }
  };

  template<class T> class vector_td<T,3>
  {
  public:

    T vec[3];

    __inline__ __host__ __device__ vector_td(){}

    __inline__ __host__ __device__ vector_td(T x, T y,T z){
      vec[0]=x;
      vec[1]=y;
      vec[2]=z;
    }

    __inline__ __host__ __device__ vector_td(T x){
      vec[0]=x;
      vec[1]=x;
      vec[2]=x;
    }

    __inline__ __host__ __device__ T& operator[](const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const long long i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const long long i) const {
      return vec[i];
    }
  };

  template<class T> class vector_td<T,4>
  {
  public:

    T vec[4];

    __inline__ __host__ __device__ vector_td(){}

    __inline__ __host__ __device__ vector_td(T x, T y,T z,T w){
      vec[0]=x;
      vec[1]=y;
      vec[2]=z;
      vec[3]=w;
    }

    __inline__ __host__ __device__ vector_td(T x){
      vec[0]=x;
      vec[1]=x;
      vec[2]=x;
      vec[3]=x;
    }

    __inline__ __host__ __device__ T& operator[](const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator[](const long long i) const {
      return vec[i];
    }

    __inline__ __host__ __device__ T& operator()(const long long i){
      return vec[i];
    }

    __inline__ __host__ __device__ const T& operator()(const long long i) const {
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

  typedef vector_td<long long,1> intd1;
  typedef vector_td<long long,2> intd2;
  typedef vector_td<long long,3> intd3;
  typedef vector_td<long long,4> intd4;

  typedef vector_td<unsigned long long,1> uintd1;
  typedef vector_td<unsigned long long,2> uintd2;
  typedef vector_td<unsigned long long,3> uintd3;
  typedef vector_td<unsigned long long,4> uintd4;

    //template <class T, unsigned long long D> std::ostream & operator<< (std::ostream & os, const vector_td<T, D>& vec)
    //{
    //    unsigned long long ii;
    //    for ( ii=0; ii<D; ii++ )
    //    {
    //        os << vec[ii] << " ";
    //    }

    //    return os;
    //}
}
