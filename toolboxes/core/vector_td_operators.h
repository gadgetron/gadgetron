/** \file vector_td_operators.h
    \brief Common operators for the vector_td class
*/

#pragma once

#include "vector_td.h"
#include "core_defines.h"

namespace Gadgetron{

  //
  // Return types
  //

  template <class T, class I> struct vectorTDReturnType {};
  template <class T> struct vectorTDReturnType<T,T> {typedef T type;};
  template<> struct vectorTDReturnType<unsigned int, int> {typedef int type;};
  template<> struct vectorTDReturnType<int, unsigned int> {typedef int type;};
  template<> struct vectorTDReturnType<int, bool> {typedef int type;};
  template<> struct vectorTDReturnType<bool,int> {typedef int type;};
  template<> struct vectorTDReturnType<unsigned int, bool> {typedef int type;};
  template<> struct vectorTDReturnType<bool,unsigned int> {typedef int type;};
  template<> struct vectorTDReturnType<float, unsigned int> {typedef float type;};
  template<> struct vectorTDReturnType<unsigned int, float> {typedef float type;};
  template<> struct vectorTDReturnType<float, int> {typedef float type;};
  template<> struct vectorTDReturnType<int, float> {typedef float type;};
  template<> struct vectorTDReturnType<float, bool> {typedef float type;};
	template<> struct vectorTDReturnType<bool, float> {typedef float type;};
  template<> struct vectorTDReturnType<double, unsigned int> {typedef double type;};
  template<> struct vectorTDReturnType<unsigned int, double> {typedef double type;};
  template<> struct vectorTDReturnType<double, int> {typedef double type;};
  template<> struct vectorTDReturnType<int, double> {typedef double type;};
  template<> struct vectorTDReturnType<double, bool> {typedef double type;};
  template<> struct vectorTDReturnType<bool, double> {typedef double type;};
  template<> struct vectorTDReturnType<double, float> {typedef double type;};
  template<> struct vectorTDReturnType<float,double> {typedef double type;};

  //
  // Operators are defined as component wise operations.
  //

  //
  // Arithmetic operators
  //

  template< class T,class R,  unsigned int D > __inline__ __host__ __device__
  void operator+= ( vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] += v2.vec[i];
  }

  template< class T,class R,  unsigned int D > __inline__ __host__ __device__
  void operator+= ( vector_td<T,D> &v1, const R &v2 )
  {
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] += v2;
  }

  template< class T,class R,  unsigned int D > __inline__ __host__ __device__
  void operator-= ( vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] -= v2.vec[i];
  }


  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  void operator*= ( vector_td<T,D> &v1, const R &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] *= v2;
  }

  template< class T,class R,  unsigned int D > __inline__ __host__ __device__
  void operator *=  ( vector_td<T,D> &v1, const vector_td<R,D> &v2 )
	{
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] *= v2.vec[i];
  }

  template< class T,class R,  unsigned int D > __inline__ __host__ __device__
  void operator /= ( vector_td<T,D> &v1, const R &v2 )
  {
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] /= v2;
  }

  template< class T,class R,  unsigned int D > __inline__ __host__ __device__
  void operator /=  ( vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  {
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] /= v2.vec[i];
  }

  template< class T,class R,  unsigned int D > __inline__ __host__ __device__
  void component_wise_div_eq ( vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] /= v2.vec[i];
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator+ ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]+v2.vec[i];
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator+ ( const vector_td<T,D> &v1, const R &v2 )
  {
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]+v2;
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator- ( const vector_td<T,D> &v1, const R &v2 )
  {
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]-v2;
    return res;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator+ (const R &v2, const vector_td<T,D> &v1 )
  {
    return v1+v2;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator- ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]-v2.vec[i];
    return res;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  vector_td<T,D> operator- ( const vector_td<T,D> &v1)
  {
    vector_td<T,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = -v1.vec[i];
    return res;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> component_wise_mul ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]*v2.vec[i];
    return res;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  vector_td<T,D> component_wise_mul ( const vector_td<T,D> &v1, const vector_td<T,D> &v2 )
  {
    vector_td<T,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]*v2.vec[i];
    return res;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator* ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  {
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ )  res.vec[i] = v1.vec[i]*v2.vec[i];
    return res;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator* ( const vector_td<T,D> &v1, const R &v2 )
  { 
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]*v2;
    return res;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator* ( const R &v1, const vector_td<T,D> &v2 )
  { 
    return v2*v1;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator/ ( const vector_td<T,D> &v1, const R &v2 )
  {
    vector_td<typename vectorTDReturnType<T,R>::type,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = v1.vec[i]/v2;
    return res;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> operator/ ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  {
    vector_td<typename vectorTDReturnType<T,R>::type,D> res = v1;
    for(unsigned int i=0; i<D; i++ ) res[i] /= v2[i];
    return res;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  vector_td<typename vectorTDReturnType<T,R>::type,D> component_wise_div ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    return v1/v2;
  }

  // 
  // "Strong" comparison operators
  //

  template< class T, unsigned int D > __inline__ __host__ __device__
  bool operator== ( const vector_td<T,D> &v1, const vector_td<T,D> &v2 ) 
  { 
    for(unsigned int i=0; i<D; i++ ) if(!(v1.vec[i] == v2.vec[i])) return false;
    return true;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  bool operator!= ( const vector_td<T,D> &v1, const vector_td<T,D> &v2 ) 
  { 
    for(unsigned int i=0; i<D; i++ ) if((v1.vec[i] != v2.vec[i])) return true;
    return false;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  bool operator&& ( const vector_td<T,D> &v1, const vector_td<T,D> &v2 ) 
  { 
    for(unsigned int i=0; i<D; i++ ) if(!(v1.vec[i] && v2.vec[i])) return false;
    return true;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  bool operator|| ( const vector_td<T,D> &v1, const vector_td<T,D> &v2 ) 
  { 
    for(unsigned int i=0; i<D; i++ ) if(!(v1.vec[i] || v2.vec[i])) return false;
    return true;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  bool operator< ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(!(v1.vec[i] < v2.vec[i])) return false;
    return true;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  bool operator<= ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(!(v1.vec[i] <= v2.vec[i])) return false;
    return true;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool operator> ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(!(v1.vec[i] > v2.vec[i])) return false;
    return true;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool operator>= ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(!(v1.vec[i] >= v2.vec[i])) return false;
    return true;
  }

  //
  // "Weak" comparison "operators"
  //

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_equal ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] == v2.vec[i]) return true;
    return false;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_not_equal ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] != v2.vec[i]) return true;
    return false;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_and ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] && v2.vec[i]) return true;
    return false;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_or ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] || v2.vec[i]) return true;
    return false;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_less ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] < v2.vec[i]) return true;
    return false;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_less_equal ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] <= v2.vec[i]) return true;
    return false;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_greater ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] > v2.vec[i]) return true;
    return false;
  }

  template< class T, class R, unsigned int D > __inline__ __host__ __device__
  bool weak_greater_equal ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    for(unsigned int i=0; i<D; i++ ) if(v1.vec[i] >= v2.vec[i]) return true;
    return false;
  }

  //
  // Vector comparison "operators"
  //

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<bool,D> vector_equal ( const vector_td<T,D> &v1, const vector_td<T,D> &v2 )
  { 
    vector_td<T,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] == v2.vec[i]);
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<bool,D> vector_not_equal ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    vector_td<T,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] != v2.vec[i]);
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<T,D> vector_and ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
    vector_td<T,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] && v2.vec[i]);
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<bool,D> vector_or ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
  	vector_td<bool,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] || v2.vec[i]);
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<bool,D> vector_less ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
  	vector_td<bool,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] < v2.vec[i]);
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<bool,D> vector_less_equal ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  { 
  	vector_td<bool,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] <= v2.vec[i]);
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<bool,D> vector_greater ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  {
    vector_td<bool,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] > v2.vec[i]);
    return res;
  }

  template< class T,class R, unsigned int D > __inline__ __host__ __device__
  vector_td<bool,D> vector_greater_equal ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
  {  
  	vector_td<bool,D> res;
    for(unsigned int i=0; i<D; i++ ) res.vec[i] = (v1.vec[i] >= v2.vec[i]);
    return res;
  }

  //
  // Integer only operators
  //

  template< class T, unsigned int D > __inline__ __host__ __device__
  void operator<<= ( vector_td<T,D> &v1, size_t shifts ) 
  { 
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] <<= shifts;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  void operator>>= ( vector_td<T,D> &v1, size_t shifts ) 
  { 
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] >>= shifts;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  vector_td<T,D> operator<< ( const vector_td<T,D> &v1, size_t shifts ) 
  { 
    vector_td<T,D> res = v1;
    res <<= shifts;
    return res;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  vector_td<T,D> operator>> ( const vector_td<T,D> &v1, size_t shifts ) 
  { 
    vector_td<T,D> res = v1;
    res >>= shifts;
    return res;
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  void operator%= ( vector_td<T,D> &v1, const vector_td<T,D> &v2 ) 
  { 
    for(unsigned int i=0; i<D; i++ ) v1.vec[i] %= v2.vec[i];
  }

  template< class T, unsigned int D > __inline__ __host__ __device__
  vector_td<T,D> operator% ( const vector_td<T,D> &v1, const vector_td<T,D> &v2 ) 
  { 
    vector_td<T,D> res = v1;
    res %= v2;
    return res;
  }


}

namespace std {

template<class T, unsigned int D> T* begin(Gadgetron::vector_td<T,D>& v){
  return v.vec;
};

template<class T, unsigned int D> T* end(Gadgetron::vector_td<T,D>& v){
return v.vec+D;
};

}
