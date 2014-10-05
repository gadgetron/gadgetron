/*
 * matrix_vector_op.h
 *
 *  Created on: Dec 9, 2011
 *      Author: hansenms
 */

#ifndef MATRIX_VECTOR_OP_H_
#define MATRIX_VECTOR_OP_H_

#include <hoNDArray.h>
#include "complext.h"

#include "linalg_export.h"

namespace Gadgetron
{

  template <typename T> EXPORTLINALG double hoNDArray_norm2(hoNDArray<T>* X);
  template <typename T> EXPORTLINALG double hoNDArray_asum(hoNDArray<T>* X);
  template <typename T> EXPORTLINALG void hoNDArray_scal(T SA, hoNDArray<T>* X);

/**
 *
 *  Performs C = alpha*(A*B) + beta*C
 *
 */
template <typename T> EXPORTLINALG void hoNDArray_gemm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha,  hoNDArray<T>* C, T beta);


/**
 *
 *  Performs Y = A*X+Y
 *
 */
 template <typename T> EXPORTLINALG void hoNDArray_axpy( T* A, hoNDArray<T>* X, hoNDArray<T>* Y); 

/**
 *  Performs B = alpha*A*B
 *
 *  A should be lower triangular.
 *
 */
template <typename T> EXPORTLINALG void hoNDArray_trmm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha);


 void elementWiseMultiply(int n,  std::complex<float> *a, std::complex<float> *x, std::complex<float> *y);


} //namespace gadgetron

#endif /* MATRIX_VECTOR_OP_H_ */
