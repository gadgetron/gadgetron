#pragma once
#include "complext.h"
#include "cublas_v2.h"
#include "ndarray_vector_td_utilities.h"


template<class T> int cublas_axpy(cublasHandle_t hndl, int n, const T* a , const T* x , int incx,  T* y, int incy);
template<class T> int cublas_dot(cublasHandle_t, int, const T*, int, const  T*, int, T*);
template<class T> int cublas_nrm2(cublasHandle_t, int, const T*, int, typename realType<T>::type *);
template<class T> int cublas_amax(cublasHandle_t handle, int n,const T *x, int incx, int *result);
template<class T> int cublas_amin(cublasHandle_t handle, int n,const T *x, int incx, int *result);
template<class T> int cublas_asum(cublasHandle_t handle, int n,const T *x, int incx, typename realType<T>::type *result);

template<class T> T dot(cuNDArray<T> *x,cuNDArray<T> *y,int device=CUNDA_NDARRAY_DEVICE );
template<class T> typename realType<T>::type nrm2( cuNDArray<T>* arr, int device=CUNDA_NDARRAY_DEVICE );

template<class T, class R> void axpy(R a, cuNDArray<T>* x, cuNDArray<T>* y, int device=CUNDA_NDARRAY_DEVICE);

template<class T> int amin(cuNDArray<T>* x,int device=CUNDA_NDARRAY_DEVICE );
template<class T> int amax(cuNDArray<T>* x,int device=CUNDA_NDARRAY_DEVICE );

template<class T> typename realType<T>::type asum(cuNDArray<T>* x,int device=CUNDA_NDARRAY_DEVICE );

