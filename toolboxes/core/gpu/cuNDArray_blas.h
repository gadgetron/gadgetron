/** \file cuNDArray_blas.h
    \brief BLAS level-1 functions on the cuNDArray class.
    
    cuNDArray_blas.h provides BLAS level-1 functions on the cuNDArray class.
    The cuNDArray is temporarily reshaped to a column vector for the respective operations.
    The implementation is based on CUBLAS.
    This code is purposely split into a header and underlying implementation (.cpp) 
    as this allows specific instantiation of the supported template types.     
    The supported types are float, double, std::complex<float>, std::complex<double>, 
    Gadgetron::complext<float>, and Gadgetron::complext<double>.
*/

#pragma once

#include "cuNDArray.h"
#include "complext.h"
#include "gpucore_export.h"

#include <cublas_v2.h>

namespace Gadgetron{

  template<class T> EXPORTGPUCORE cublasStatus_t cublas_axpy(cublasHandle_t hndl, int n, const T* a , const T* x , int incx,  T* y, int incy);
  template<class T> EXPORTGPUCORE cublasStatus_t cublas_dot(cublasHandle_t, int, const T*, int, const  T*, int, T*, bool cc = true);
  template<class T> EXPORTGPUCORE cublasStatus_t cublas_nrm2(cublasHandle_t, int, const T*, int, typename realType<T>::Type *result);
  template<class T> EXPORTGPUCORE cublasStatus_t cublas_amax(cublasHandle_t handle, int n,const T *x, int incx, int *result);
  template<class T> EXPORTGPUCORE cublasStatus_t cublas_amin(cublasHandle_t handle, int n,const T *x, int incx, int *result);
  template<class T> EXPORTGPUCORE cublasStatus_t cublas_asum(cublasHandle_t handle, int n,const T *x, int incx, typename realType<T>::Type *result);

  template<class T> EXPORTGPUCORE T dot( cuNDArray<T> *x, cuNDArray<T> *y, bool cc = true );

  template<class T> EXPORTGPUCORE typename realType<T>::Type nrm2( cuNDArray<T> *x );

  template<class T> EXPORTGPUCORE void axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y );

  template<class T> EXPORTGPUCORE void axpy( T a, cuNDArray<complext<T> > *x, cuNDArray<complext<T> > *y );
  
  /**
   * @brief Gets the index of the index of the element with minimum absolute
   * @param x Input data
   * @return index of absolute minimum values
   */
  template<class T> EXPORTGPUCORE int amin( cuNDArray<T> *x );
  
  /**
   * @brief Gets the index of the index of the element with maximum absolute
   * @param x Input data
   * @return index of absolute maximum values
   * @details Note that this returns the C-style index and NOT the Fortran index.
   */
  template<class T> EXPORTGPUCORE int amax( cuNDArray<T> *x);
  
  template<class T> EXPORTGPUCORE typename realType<T>::Type asum( cuNDArray<T> *x );
  
  EXPORTGPUCORE std::string getCublasErrorString(cublasStatus_t err);
}
