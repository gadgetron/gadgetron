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

//Declaration of BLAS routines
/*
 * We will opt to not use the easier CBLAS interface to give us the best change of being compatible on all platforms.
 * We will declare the BLAS (and LAPACK) routines ourselves.
 *
 */
extern "C" {
	//GEMM - Generalized matrix-matrix multiplication
	void sgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);
	void dgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);
	void cgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);
	void zgemm_(char* TRANSA,char* TRANSB,int* M, int *N, int *K, void* ALPHA,
				void *A, int* LDA, void* B, int* LDB, void* BETA, void* C, int *LDC);

	//TRMM - Multiplication with a triangular matrix
	void strmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M,int* N,
			void* ALPHA,void* A,int* LDA,void* B, int* LDB);
	void dtrmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M,int* N,
			void* ALPHA,void* A,int* LDA,void* B, int* LDB);
	void ctrmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M,int* N,
			void* ALPHA,void* A,int* LDA,void* B, int* LDB);
	void ztrmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M,int* N,
			void* ALPHA,void* A,int* LDA,void* B, int* LDB);
}


namespace Gadgetron
{

/**
 *
 *  Performs C = alpha*(A*B) + beta*C
 *
 */
template <typename T> EXPORTLINALG void hoNDArray_gemm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha,  hoNDArray<T>* C, T beta);

/**
 *  Performs B = alpha*A*B
 *
 *  A should be lower triangular.
 *
 */
template <typename T> EXPORTLINALG void hoNDArray_trmm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha);

} //namespace gadgetron

#endif /* MATRIX_VECTOR_OP_H_ */
