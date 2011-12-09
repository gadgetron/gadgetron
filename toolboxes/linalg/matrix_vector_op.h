/*
 * matrix_vector_op.h
 *
 *  Created on: Dec 9, 2011
 *      Author: hansenms
 */

#ifndef MATRIX_VECTOR_OP_H_
#define MATRIX_VECTOR_OP_H_

#include <hoNDArray.h>
#include <complex>

#include "linalg_export.h"

#if defined __APPLE__
	#include <Accelerate/Accelerate.h>
#else
	#include <cblas.h>
#endif


/**
 *
 *  Performs C = alpha*(A*B) + beta*C
 *
 */
template <typename T> EXPORTLINALG int hoNDArray_gemm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha,  hoNDArray<T>* C, T beta);


#endif /* MATRIX_VECTOR_OP_H_ */
