#pragma once

#include "vector_td.h"
#include "vector_td_operators.h"

#include <float.h>
#include <cmath>

//
// Get/set operations on vector_td<T,D>
//

template<class T, unsigned int D> __inline__ __host__ __device__ T get( const vector_td<T,D> vec, unsigned int dim ) { return vec.vec[dim]; }
template<class T, unsigned int D> __inline__ __host__ __device__ void set( vector_td<T,D> &vec, unsigned int dim, T val ) { vec.vec[dim] = val; }

//
// Operations on intergers/float/double and complext<REAL>
//

template<class T> __inline__ __host__ __device__ T get_zero();
template<class T> __inline__ __host__ __device__ T get_one();

//
// Operations on float/double and complext<REAL>
//

template<class T> __inline__ __host__ __device__ T get_half();

template<class T> __inline__ __host__ __device__ T reciprocal( const T );

template<class T> __inline__ __host__ __device__ T mul( const T, const T );
template<class REAL> __inline__ __host__ __device__ typename complext<REAL>::Type mul( const REAL, const typename complext<REAL>::Type );
template<class REAL> __inline__ __host__ __device__ typename complext<REAL>::Type mul( const typename complext<REAL>::Type, const REAL );

//
// Component-wise math operations
//

template<class T, unsigned int D> __inline__ __host__ __device__ vector_td<T,D> abs( const vector_td<T,D> vec )
{
  vector_td<T,D> res;
  for (unsigned int i=0; i<D; i++) {
    res.vec[i] = abs(vec.vec[i]);
  }
  return res;
}

template<class REAL, unsigned int D> __inline__ __host__ __device__ vector_td<REAL,D> ceil( const vector_td<REAL,D> vec )
{
  vector_td<REAL,D> res;
  for (unsigned int i=0; i<D; i++) {
    res.vec[i] = ceil(vec.vec[i]);
  }
  return res;
}

template<class REAL, unsigned int D> __inline__ __host__ __device__ vector_td<REAL,D> floor( const vector_td<REAL,D> vec )
{
  vector_td<REAL,D> res;
  for (unsigned int i=0; i<D; i++) {
    res.vec[i] = floor(vec.vec[i]);
  }
  return res;
}

//
// Vectorize a scalar value
//

template<class T, unsigned int D> __inline__ __host__ __device__ vector_td<T,D> to_vector_td( const T scalar )
{
  vector_td<T,D> res;
  for (unsigned int i=0; i<D; i++) {
    res.vec[i] = scalar;
  }
  return res;
}

//
// Grid <-> index transformations
//

template<unsigned int D> __inline__ __host__ __device__ typename uintd<D>::Type idx_to_co( unsigned int idx, const typename uintd<D>::Type dims )
{
  typename uintd<D>::Type co;
  unsigned int idx_tmp = idx;
  for (unsigned int i=0; i<D; i++) {
    co.vec[i] = idx_tmp%dims.vec[i];
    idx_tmp -= co.vec[i];
    idx_tmp /= dims.vec[i];
  }
  return co;
} 

template<unsigned int D> __inline__ __host__ __device__ typename intd<D>::Type idx_to_co( int idx, const typename intd<D>::Type dims )
{
  typename intd<D>::Type co;
  int idx_tmp = idx;
  for (unsigned int i=0; i<D; i++) {
    co.vec[i] = idx_tmp%dims.vec[i];
    idx_tmp -= co.vec[i];
    idx_tmp /= dims.vec[i];
  }
  return co;
} 

template<unsigned int D> __inline__ __host__ __device__ unsigned int co_to_idx( const typename uintd<D>::Type co, const typename uintd<D>::Type dims )
{
  unsigned int idx = 0;
  unsigned long block_size = 1;
  for (unsigned int i=0; i<D; i++) {
    idx += (block_size*co.vec[i]);
    block_size *= dims.vec[i];
  }
  return idx;
} 

template<unsigned int D> __inline__ __host__ __device__ unsigned int co_to_idx( const typename uintd<D>::Type co, const typename uintd<D>::Type dims, const typename uintd<D>::Type &order )
{
  unsigned int idx = 0;
  unsigned long block_size = 1;
  for (unsigned int i=0; i<D; i++){
    idx += (block_size*co.d[order.vec[i]]);
    block_size *= dims.d[order.vec[i]];
  }
  return idx;
} 

template<unsigned int D> __inline__ __host__ __device__ int co_to_idx( const typename intd<D>::Type co, const typename intd<D>::Type dims )
{
  int idx = 0;
  unsigned long block_size = 1;
  for (unsigned int i=0; i<D; i++) {
    idx += (block_size*co.vec[i]);
    block_size *= dims.vec[i];
  }
  return idx;
} 

template<unsigned int D> __inline__ __host__ __device__ int co_to_idx( const typename intd<D>::Type co, const typename intd<D>::Type dims, const typename intd<D>::Type order )
{
  int idx = 0;
  unsigned long block_size = 1;
  for (unsigned int i=0; i<D; i++){
    idx += (block_size*co.d[order.vec[i]]);
    block_size *= dims.d[order.vec[i]];
  }
  return idx;
} 

template<unsigned int D> __inline__ __host__ __device__ typename uintd<D>::Type counting_vec()
{
  typename uintd<D>::Type res;
  for(unsigned int i=0; i<D; i++) {
    res.vec[i]=i;
  }
  return res;
}

//
// Conversion between vector<unsigned int> and uintd
//

template<unsigned int D> 
std::vector<unsigned int> uintd_to_vector( typename uintd<D>::Type _uintd )
{
  std::vector<unsigned int> out(D);
  for( unsigned int i=0; i<D; i++ )
    out[i] = _uintd.vec[i];
  return out;
}

template<unsigned int D> 
typename uintd<D>::Type vector_to_uintd( std::vector<unsigned int> _vector )
{
  typename uintd<D>::Type out;
  if( _vector.size() < D ){
    std::cout << "Cannot convert vector to typename uintd<D>" << std::endl;
    return typename uintd<D>::Type();
  }
  
  for( unsigned int i=0; i<D; i++ ){
    out.vec[i] = _vector[i];
  }
  
  return out;
}

//
// Reductions on vector_td<T,D>
//

template<class T, unsigned int D> __inline__ __host__ __device__ T prod( const vector_td<T,D> vec ){
  T res = get_one<T>();
  for (unsigned int i=0; i<D; i++){
    res *= vec.vec[i];
  }
  return res;
}

template<class T, unsigned int D> __inline__ __host__ __device__ T sum( const vector_td<T,D> vec ){
  T res = vec.vec[0];
  for (unsigned int i=1; i<D; i++){
    res += vec.vec[i];
  }
  return res;
}

template<class T, unsigned int D> __inline__ __host__ __device__ T dot( const vector_td<T,D> vec1, const vector_td<T,D> vec2 ){
  T res = (vec1.vec[0]*vec2.vec[0]);
  for (unsigned int i=1; i<D; i++){
    res += (vec1.vec[i]*vec2.vec[i]);
  }
  return res;
}

//
// Reductions on reald<REAL,D>
//

template<class REAL, unsigned int D> __inline__ __host__ __device__ REAL norm_squared( const typename reald<REAL,D>::Type vec ){
  REAL res = get_zero<REAL>();
  for (unsigned int i=0; i<D; i++){
    res += (vec.vec[i]*vec.vec[i]);
  }

  return res;
}

template<class REAL, unsigned int D> __inline__ __host__ __device__ REAL norm( const typename reald<REAL,D>::Type vec ){
  return sqrt(norm_squared<REAL,D>(vec));
}

//
// Reductions on float/double and complext<REAL>
//

template<class REAL, class T> __inline__ __host__ __device__ REAL norm_squared( const T );

template<class REAL, class T> __inline__ __host__ __device__ REAL norm( const T );


//
// Type conversion
//

template<class T, unsigned int D> __inline__ __host__ __device__ vector_td<int,D> to_intd( const vector_td<T,D> vec )
{
  vector_td<int,D> res;
  for (unsigned int i=0; i<D; i++){
    res.vec[i] = (int) vec.vec[i];
  }
  return res;
}

template<class T, unsigned int D> __inline__ __host__ __device__ vector_td<unsigned int,D> to_uintd( const vector_td<T,D> vec )
{
  vector_td<unsigned int,D> res;
  for (unsigned int i=0; i<D; i++){
    res.vec[i] = (unsigned int) vec.vec[i];
  }
  return res;
}

template<class REAL, class T, unsigned int D> __inline__ __host__ __device__ vector_td<REAL,D> to_reald( const vector_td<T,D> vec )
{
  vector_td<REAL,D> res;
  for (unsigned int i=0; i<D; i++){
    res.vec[i] = (REAL) vec.vec[i];
  }
  return res;
}

//
// Trigonometry on complext<REAL>
//

template<class REAL, class T> __inline__ __host__ __device__ REAL real( const T );

template<class REAL> __inline__ __host__ __device__ REAL imag( const typename complext<REAL>::Type z ){
  return z.vec[1];
}

template<class REAL> __inline__ __host__ __device__ REAL arg( const typename complext<REAL>::Type z ){
  return atan2(imag<REAL>(z), real<REAL>(z));
}

template<class T> __inline__ __host__ __device__ T conj( const T );

//
// Explicit instantiations 
//

template<> __inline__ __host__ __device__ float get_zero<float>()
{
  return 0.0f;
}

template<> __inline__ __host__ __device__ double get_zero<double>()
{
  return 0.0;
}

template<> __inline__ __host__ __device__ int get_zero<int>()
{
  return 0;
}

template<> __inline__ __host__ __device__ unsigned int get_zero<unsigned int>()
{
  return 0;
}

template<> __inline__ __host__ __device__ vector_td<float,2> get_zero<vector_td<float,2> >()
{
  vector_td<float,2> res;
  res.vec[0] = 0.0f;
  res.vec[1] = 0.0f;
  return res;
}

template<> __inline__ __host__ __device__ vector_td<double,2> get_zero<vector_td<double,2> >()
{
  vector_td<double,2> res;
  res.vec[0] = 0.0;
  res.vec[1] = 0.0;
  return res;
}

template<> __inline__ __host__ __device__ float get_one<float>()
{
  return 1.0f;
}

template<> __inline__ __host__ __device__ double get_one<double>()
{
  return 1.0;
}

template<> __inline__ __host__ __device__ int get_one<int>()
{
  return 1;
}

template<> __inline__ __host__ __device__ unsigned int get_one<unsigned int>()
{
  return 1;
}

template<> __inline__ __host__ __device__ vector_td<float,2> get_one<vector_td<float,2> >()
{
  vector_td<float,2> res;
  res.vec[0] = 1.0f;
  res.vec[1] = 0.0f;
  return res;
}

template<> __inline__ __host__ __device__ vector_td<double,2> get_one<vector_td<double,2> >()
{
  vector_td<double,2> res;
  res.vec[0] = 1.0;
  res.vec[1] = 0.0;
  return res;
}

template<> __inline__ __host__ __device__ float get_half<float>()
{
  return 0.5f;
}

template<> __inline__ __host__ __device__ double get_half<double>()
{
  return 0.5;
}

template<> __inline__ __host__ __device__ vector_td<float,2> get_half<vector_td<float,2> >()
{
  vector_td<float,2> res;
  res.vec[0] = 0.5f;
  res.vec[1] = 0.0f;
  return res;
}

template<> __inline__ __host__ __device__ vector_td<double,2> get_half<vector_td<double,2> >()
{
  vector_td<double,2> res;
  res.vec[0] = 0.5;
  res.vec[1] = 0.0;
  return res;
}

template <> __inline__ __host__ __device__ float norm_squared<float,float>( const float v ){
  return v*v;
}

template <> __inline__ __host__ __device__ double norm_squared<double,double>( const double v ){
  return v*v;
}

template <> __inline__ __host__ __device__ float norm_squared<float,float_complext::Type>( const float_complext::Type z ){
  return norm_squared<float,2>(z);
}

template <> __inline__ __host__ __device__ double norm_squared<double,double_complext::Type>( const double_complext::Type z ){
  return norm_squared<double,2>(z);
}

template<> __inline__ __host__ __device__ float norm<float,float>( const float v ){
  return abs(v);
}

template<> __inline__ __host__ __device__ double norm<double,double>( const double v ){
  return abs(v);
}

template<> __inline__ __host__ __device__ float norm<float,float_complext::Type>( const float_complext::Type z ){
  return norm<float,2>(z);
}

template<> __inline__ __host__ __device__ double norm<double,double_complext::Type>( const double_complext::Type z ){
  return norm<double,2>(z);
}

template<> __inline__ __host__ __device__ float reciprocal<float>( const float real )
{
  return 1.0f/real;
}

template<> __inline__ __host__ __device__ double reciprocal<double>( const double real )
{
  return 1.0/real;
}

template<> __inline__ __host__ __device__ float conj<float>( const float r ){
  return r;
}

template<> __inline__ __host__ __device__ double conj<double>( const double r ){
  return r;
}

template<> __inline__ __host__ __device__ float_complext::Type conj<float_complext::Type>( const float_complext::Type z ){
  float_complext::Type res;
  res.vec[0] = z.vec[0]; 
  res.vec[1] = -z.vec[1];
  return res;
}

template<> __inline__ __host__ __device__ double_complext::Type conj<double_complext::Type>( const double_complext::Type z ){
  double_complext::Type res;
  res.vec[0] = z.vec[0]; 
  res.vec[1] = -z.vec[1];
  return res;
}

template<> __inline__ __host__ __device__ float_complext::Type reciprocal<float_complext::Type>( const float_complext::Type z )
{
  float_complext::Type res = conj<float_complext::Type>(z);
  res *= (1.0f/norm_squared<float,float_complext::Type>(z));
  return res;
}

template<> __inline__ __host__ __device__ double_complext::Type reciprocal<double_complext::Type>( const double_complext::Type z )
{
  double_complext::Type res = conj<double_complext::Type>(z);
  res *= (1.0/norm_squared<double,double_complext::Type>(z));
  return res;
}

template<> __inline__ __host__ __device__ float mul<float>( const float a, const float b )
{
  return a*b;
}

template<> __inline__ __host__ __device__ float_complext::Type mul<float_complext::Type>( const float_complext::Type a, const float_complext::Type b )
{
  float_complext::Type res;
  res.vec[0] = a.vec[0]*b.vec[0]-a.vec[1]*b.vec[1];
  res.vec[1] = a.vec[0]*b.vec[1]+a.vec[1]*b.vec[0];
  return res;
}

template<> __inline__ __host__ __device__ float_complext::Type mul<float>( const float a, const float_complext::Type b )
{
  float_complext::Type res;
  res.vec[0] = a*b.vec[0];
  res.vec[1] = a*b.vec[1];
  return res;
}

template<> __inline__ __host__ __device__ float_complext::Type mul<float>( const float_complext::Type a, const float b )
{
  float_complext::Type res;
  res.vec[0] = b*a.vec[0];
  res.vec[1] = b*a.vec[1];
  return res;
}

template<> __inline__ __host__ __device__ double mul<double>( const double a, const double b )
{
  return a*b;
}

template<> __inline__ __host__ __device__ double_complext::Type mul<double_complext::Type>( const double_complext::Type a, const double_complext::Type b )
{
  double_complext::Type res;
  res.vec[0] = a.vec[0]*b.vec[0]-a.vec[1]*b.vec[1];
  res.vec[1] = a.vec[0]*b.vec[1]+a.vec[1]*b.vec[0];
  return res;
}

template<> __inline__ __host__ __device__ double_complext::Type mul<double>( const double a, const double_complext::Type b )
{
  double_complext::Type res;
  res.vec[0] = a*b.vec[0];
  res.vec[1] = a*b.vec[1];
  return res;
}

template<> __inline__ __host__ __device__ double_complext::Type mul<double>( const double_complext::Type a, const double b )
{
  double_complext::Type res;
  res.vec[0] = b*a.vec[0];
  res.vec[1] = b*a.vec[1];
  return res;
}

template<> __inline__ __host__ __device__ float real<float,float>( float r ){
  return r;
}

template<> __inline__ __host__ __device__ double real<double,double>( double r ){
  return r;
}

template<> __inline__ __host__ __device__ float real<float,float_complext::Type>( const float_complext::Type z ){
  return z.vec[0];
}

template<> __inline__ __host__ __device__ double real<double,double_complext::Type>( const double_complext::Type z ){
  return z.vec[0];
}
