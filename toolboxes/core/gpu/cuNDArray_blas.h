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


#include <cublas_v2.h>

namespace Gadgetron{

    /**
   * @brief compute dot product of two vectors
   * @param x Input data
   * @param y Input data
   * @param batchSize Size of each batch of the data that is processed on the GPU at a time 
   * @param cc complex conjugate. calculate dot product of conj(x) and y if cc = true else x.y
   * @return x.y
   */

  template<class T> T dot( cuNDArray<T> *x, cuNDArray<T> *y, size_t batchSize = INT_MAX, bool cc = true );

  
  /**
   * @brief compute ||x||
   * @param x Input data
   * @param batchSize Size of each batch of the data that is processed on the GPU at a time 
   * @return Euclidean norm of x
   */

  template<class T> typename realType<T>::Type nrm2( cuNDArray<T> *x , size_t batchSize = INT_MAX);

    /**
   * @brief Compute  y = ax + y
   * @param a scalar a
   * @param x Input data
   * @param y Input data that will be updated
   * @param batchSize Size of each batch of the data that is processed on the GPU at a time
   */

  template<class T> void axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y , size_t batchSize = INT_MAX);

  template<class T> void axpy( T a, cuNDArray<complext<T> > *x, cuNDArray<complext<T> > *y , size_t batchSize = INT_MAX);
  
  /**
   * @brief Gets the index of the index of the element with minimum absolute
   * @param x Input data
   * @param batchSize Size of each batch of the data that is processed on the GPU at a time 
   * @return index of absolute minimum values
   * @details Note that this returns the C-style index and NOT the Fortran index.
   */

  template<class T> size_t amin( cuNDArray<T> *x , size_t batchSize = INT_MAX);

  /**
   * @brief Gets the index of the index of the element with maximum absolute
   * @param x Input data
   * @param batchSize Size of each batch of the data that is processed on the GPU at a time
   * @return index of absolute maximum values
   * @details Note that this returns the C-style index and NOT the Fortran index.
   */

  template<class T> size_t amax( cuNDArray<T> *x , size_t batchSize = INT_MAX);

  /**
   * @brief Compute the sum of absolute values of the elements of vector x
   * @param x Input data
   * @param batchSize Size of each batch of the data that is processed on the GPU at a time
   * @return absolute value sum of the input vector
   */
  
  template<class T> typename realType<T>::Type asum( cuNDArray<T> *x , size_t batchSize = INT_MAX );
  
  std::string gadgetron_getCublasErrorString(cublasStatus_t err);

}
