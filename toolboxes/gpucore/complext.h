/*
 * complext.h
 *
 *  Created on: Feb 13, 2012
 *      Author: David C Hansen
 */

#pragma once

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
};

template<class T>  __inline__ __host__ __device__ complext<T> polar(const T& rho, const T& theta = 0){
  return complext<T>(rho*cos(theta),rho*sin(theta));
}

template<class T> __inline__ __host__ __device__ T abs(complext<T> comp){
  return sqrt(comp.vec[0]*comp.vec[0]+comp.vec[1]*comp.vec[1]);
}

template<class T> __inline__ __host__ __device__ complext<T> sin(complext<T> comp){
  return complext<T>(sin(comp.vec[0])*cosh(comp.vec[1]),cos(comp.vec[0])*sinh(comp.vec[1]));
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
  return atan2(comp.vec[1],comp.vec[0]);
}

template<class T> __inline__ __host__ __device__  complext<T> operator*(const T& r,const complext<T>& z){
  return complext<T>(z.vec[0]*r,z.vec[1]*r);
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

__inline__ __host__ __device__ double conj(const double& r) { return r;}

__inline__ __host__ __device__ float conj(const float& r) { return r;}


template<class T> __inline__ __host__ __device__ complext<T> conj( const complext<T>& z ){
  complext<T> res=z;
  res.conj();
  return res;
}

typedef complext<float> float_complext;
typedef complext<double> double_complext;
