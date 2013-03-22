/** \file complext.h
    \brief An implementation of complex numbers that works for both the cpu and gpu.

    complext.h provides an implementation of complex numbers that, unlike std::complex,
    works on both the cpu and gpu. 
    It follows the interface defined for std::complex.
*/

#pragma once

#include "core_defines.h"
#include <complex>
#include <cmath>

namespace Gadgetron{

  using std::abs; // workaround for nvcc

/** 
 * \class complext
 * \brief An implementation of complex numbers that works for both the cpu and gpu.
 */
  template< class T > class complext
  {
  public:

    T vec[2];

    __inline__ __host__ __device__  T real(){
      return vec[0];
    }

    __inline__ __host__ __device__  T imag(){
      return vec[1];
    }

    __inline__ __host__ __device__  complext() {}

    __inline__ __host__ __device__  complext(T real, T imag){
      vec[0]=real;
      vec[1]=imag;
    }

    __inline__ __host__ __device__  complext(const complext<T>& tmp){
      vec[0] = tmp.vec[0];
      vec[1] = tmp.vec[1];
    }

    __inline__ __host__ __device__  complext(const T r){
      vec[0] = r;
      vec[1] = 0;
    }

    __inline__ __host__ __device__ void conj(){
      vec[1] = -vec[1];
    }

    __inline__ __host__ __device__  complext<T> operator+(const complext<T>& other){
      return complext<T>(vec[0]+other.vec[0],vec[1]+other.vec[1]);
    }

    __inline__ __host__ __device__  complext<T> operator-(const complext<T>& other){
      return complext<T>(vec[0]-other.vec[0],vec[1]-other.vec[1]);
    }

    __inline__ __host__ __device__  complext<T> operator-(){
      return complext<T>(-vec[0],-vec[1]);
    }

    __inline__ __host__ __device__  void operator-=(const complext<T>& other){
      vec[0] -= other.vec[0];
      vec[1] -= other.vec[1];
    }

    __inline__ __host__ __device__  void operator+=(const complext<T>& other){
      vec[0] += other.vec[0];
      vec[1] += other.vec[1];
    }

    __inline__ __host__ __device__  complext<T> operator*(const T& other){
      return complext<T>(vec[0]*other,vec[1]*other);
    }

    __inline__ __host__ __device__  complext<T> operator*(const complext<T>& other){
      return complext<T>(vec[0]*other.vec[0]-vec[1]*other.vec[1],vec[0]*other.vec[1]+vec[1]*other.vec[0]);
    }

    __inline__ __host__ __device__  complext<T> operator/(const T& other){
      return complext<T>(vec[0]/other,vec[1]/other);
    }

    __inline__ __host__ __device__  complext<T> operator/(const complext<T>& other){
      T cd = other.vec[0]*other.vec[0]+other.vec[1]*other.vec[1];
      return complext<T>((vec[0]*other.vec[0]+vec[1]*other.vec[1])/cd ,(vec[1]*other.vec[0]-vec[0]*other.vec[1])/cd);
    }

    __inline__ __host__ __device__  void operator*=(const T& other){
      vec[0] *= other;
      vec[1] *= other;
    }

    __inline__ __host__ __device__  void operator*=(const complext<T>& other){
      complext<T> tmp = *this;
      vec[0] = tmp.vec[0]*other.vec[0]-tmp.vec[1]*other.vec[1];
      vec[1] = tmp.vec[0]*other.vec[1]+tmp.vec[1]*other.vec[0];
    }

    __inline__ __host__ __device__  void operator/=(const T& other){
      vec[0] /= other;
      vec[1] /= other;
    }

    __inline__ __host__ __device__  void operator/=(const complext<T>& other){
      complext<T> tmp = (*this)/other;
      vec[0]=tmp.vec[0];
      vec[1]=tmp.vec[1];
    }

    __inline__ __host__ __device__  bool operator==(const complext<T>& comp2){

      return vec[0]==comp2.vec[0] && vec[1]==comp2.vec[1];
    }
    __inline__ __host__ __device__  bool operator!=(const complext<T>& comp2){

      return not(*this==comp2);
    }
  };
  
  typedef complext<float> float_complext;
  typedef complext<double> double_complext;

  template <class T> struct realType {};
  template<> struct realType<float_complext> {typedef float Type; };
  template<> struct realType<double_complext> {typedef double Type; };
  template<> struct realType<float> {typedef float Type; };
  template<> struct realType<double> {typedef double Type; };
  template<> struct realType<std::complex<float> > {typedef float Type; };
  template<> struct realType<std::complex<double> > {typedef double Type; };

  template<class T> struct stdType {typedef T type;};
  template<> struct stdType<double_complext> {typedef std::complex<double> Type;};
  template<> struct stdType<float_complext> {typedef std::complex<float> Type;};
  template<> struct stdType<std::complex<double> > {typedef std::complex<double> Type;};
  template<> struct stdType<std::complex<float> > {typedef std::complex<float> Type;};
  template<> struct stdType<double> {typedef double Type;};
  template<> struct stdType<float> {typedef float Type;};

  template<class T>  __inline__ __host__ __device__ complext<T> polar(const T& rho, const T& theta = 0){
    return complext<T>(rho*std::cos(theta),rho*std::sin(theta));
  }

  template<class T>  __inline__ __host__ __device__ complext<T> sqrt(complext<T> x){
    T r = abs(x);
    return complext<T>(::sqrt((r+x.real())/2),sgn(x.imag())*::sqrt((r-x.real())/2));
  }

  template<class T> __inline__ __host__ __device__ T abs(complext<T> comp){
    return ::sqrt(comp.vec[0]*comp.vec[0]+comp.vec[1]*comp.vec[1]);
  }

  template<class T> __inline__ __host__ __device__ complext<T> sin(complext<T> comp){
    return complext<T>(sin(comp.vec[0])*std::cosh(comp.vec[1]),std::cos(comp.vec[0])*std::sinh(comp.vec[1]));
  }

  template<class T> __inline__ __host__ __device__ complext<T> cos(complext<T> comp){
    return complext<T>(cos(comp.vec[0])*cosh(comp.vec[1]),-sin(comp.vec[0])*sinh(comp.vec[1]));
  }

  template<class T> __inline__ __host__ __device__ T imag(complext<T> comp){
    return comp.vec[1];
  }

  __inline__ __host__ __device__ double real(double r){
    return r;
  }

  __inline__ __host__ __device__ double imag(double r){
    return 0.0;
  }

  __inline__ __host__ __device__ float real(float r){
    return r;
  }

  __inline__ __host__ __device__ float imag(float r){
    return 0.0f;
  }

  template<class T> __inline__ __host__ __device__ T real(complext<T> comp){
    return comp.vec[0];
  }

  template<class T> __inline__ __host__ __device__ T arg(complext<T> comp){
    return std::atan2(comp.vec[1],comp.vec[0]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator*(const T& r,const complext<T>& z){
    return complext<T>(z.vec[0]*r,z.vec[1]*r);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator*(const complext<T>& z,const T& r){
    return complext<T>(z.vec[0]*r,z.vec[1]*r);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator+(const complext<T>& z1,const complext<T>& z2){
    return complext<T>(z1.vec[0]+z2.vec[0],z1.vec[1]+z2.vec[1]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator+(const complext<T>& z1,const T& r){
    return complext<T>(z1.vec[0]+r, z1.vec[1]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator+(const T& r,const complext<T>& z1){
    return complext<T>(z1.vec[0]+r, z1.vec[1]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator-(const complext<T>& z1,const complext<T>& z2){
    return complext<T>(z1.vec[0]-z2.vec[0],z1.vec[1]-z2.vec[1]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator-(const T& r,const complext<T>& z2){
    return complext<T>(r-z2.vec[0],-z2.vec[1]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator-(const complext<T>& z2,const T& r){
    return complext<T>(z2.vec[0]-r,z2.vec[1]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator*(const complext<T>& z1,const complext<T>& z2){
    return complext<T>(z1.vec[0]*z2.vec[0]-z1.vec[1]*z2.vec[1],z1.vec[0]*z2.vec[1]+z1.vec[1]*z2.vec[0]);
  }

  template<class T> __inline__ __host__ __device__  complext<T> operator/(const complext<T>& z1,const complext<T>& z2){
    T cd = z2.vec[0]*z2.vec[0]+z2.vec[1]*z2.vec[1];
    return complext<T>((z1.vec[0]*z2.vec[0]+z1.vec[1]*z2.vec[1])/cd ,(z1.vec[1]*z2.vec[0]-z1.vec[0]*z2.vec[1])/cd);
  }

  template<class REAL, class T> __inline__ __host__ __device__  complext<T> operator/(const REAL& real, const complext<T>& comp){
    T cd = comp.vec[0]*comp.vec[0]+comp.vec[1]*comp.vec[1];
    return complext<T>(comp.vec[0]*real/cd,-real*comp.vec[1]/cd);
  }

  template<class REAL, class T> __inline__ __host__ __device__  complext<T> operator/(const complext<T>& comp,const REAL& real){
    return complext<T>(comp.vec[0]/real,comp.vec[1]/real);
  }

  __inline__ __host__ __device__ float norm(const float& r){
    return r*r;
  }

  __inline__ __host__ __device__ double norm(const double& r){
    return r*r;
  }

  template<class T> __inline__ __host__ __device__ T norm(const complext<T>& z){
    return z.vec[0]*z.vec[0]+z.vec[1]*z.vec[1];
  }

  __inline__ __host__ __device__ double conj(const double& r){ 
    return r; }

  __inline__ __host__ __device__ float conj(const float& r) { 
    return r; }
  
  template<class T> __inline__ __host__ __device__ complext<T> conj( const complext<T>& z ){
    complext<T> res=z;
    res.conj();
    return res;
  }

  __inline__ __host__ __device__ double sgn(double x){
    return (0 < x) - (x < 0);
  }
  __inline__ __host__ __device__ float sgn(float x){
    return (0 < x) - (x < 0);
  }

  template<class T> __inline__ __host__ __device__ complext<T> sgn(complext<T> x){
    if (norm(x) <= 0) return 0;
    return (x/abs(x));
  }
}
