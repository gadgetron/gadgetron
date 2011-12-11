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

//Declaration of BLAS routines
/*
 * We will opt to not use the easier CBLAS interface to give us the best change of being compatible on all platforms.
 * We will declare the BLAS (and LAPACK) routines ourselves.
 *
 */
extern "C" {
	//GEMM
	void sgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);
	void dgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);
	void cgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);
	void zgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);
}

/**
 *
 *  Performs C = alpha*(A*B) + beta*C
 *
 */
template <typename T> EXPORTLINALG int hoNDArray_gemm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha,  hoNDArray<T>* C, T beta);


#endif /* MATRIX_VECTOR_OP_H_ */
