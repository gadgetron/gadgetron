#pragma once
#include "complext.h"
#include "cublas_v2.h"
#include "ndarray_vector_td_utilities.h"


template<class T> cublasStatus_t cublas_axpy(cublasHandle_t hndl, int n, const T* a , const T* x , int incx,  T* y, int incy);
template<class T> cublasStatus_t cublas_dot(cublasHandle_t, int, const T*, int, const  T*, int, T*);
template<class T> cublasStatus_t cublas_nrm2(cublasHandle_t, int, const T*, int, typename realType<T>::type *);
template<class T> cublasStatus_t cublas_amax(cublasHandle_t handle, int n,const T *x, int incx, int *result);
template<class T> cublasStatus_t cublas_amin(cublasHandle_t handle, int n,const T *x, int incx, int *result);
template<class T> cublasStatus_t cublas_asum(cublasHandle_t handle, int n,const T *x, int incx, typename realType<T>::type *result);

template<class T> T dot(cuNDArray<T> *x,cuNDArray<T> *y,int device );
template<class T> T dot(cuNDArray<T> *x,cuNDArray<T> *y){ return dot(x,y,x->get_device()); }

template<class T> typename realType<T>::type nrm2( cuNDArray<T>* arr, int device);
template<class T> typename realType<T>::type nrm2( cuNDArray<T>* arr){return nrm2(arr,arr->get_device());}

template<class T> void axpy(T a, cuNDArray<T>* x, cuNDArray<T>* y, int device);
template<class T> void axpy(T a, cuNDArray<T>* x, cuNDArray<T>* y){axpy(a,x,y,x->get_device());}

template<class T> void axpy(T a, cuNDArray<complext<T> >* x, cuNDArray<complext<T> >* y, int device){
	axpy(complext<T>(a),x,y,device);
}
template<class T> void axpy(T a, cuNDArray<complext<T> >* x, cuNDArray<complext<T> >* y){axpy(a,x,y,x->get_device());}


template<class T> int amin(cuNDArray<T>* x,int device);
template<class T> int amin(cuNDArray<T>* x){return amin(x,x->get_device());}

template<class T> int amax(cuNDArray<T>* x,int device);
template<class T> int amax(cuNDArray<T>* x){return amax(x,x->get_device());}

template<class T> typename realType<T>::type asum(cuNDArray<T>* x,int device);
template<class T> typename realType<T>::type asum(cuNDArray<T>* x){return asum(x,x->get_device());}

std::string getCublasErrorString(cublasStatus_t err);

