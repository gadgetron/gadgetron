/*
 * matrix_decomposition.h
 *
 *  Created on: Dec 10, 2011
 *      Author: Michael S. Hansen
 */

#ifndef MATRIX_DECOMPOSITION_H_
#define MATRIX_DECOMPOSITION_H_

#include "hoNDArray.h"

#include "linalg_export.h"

//Declaration of lapack routines
extern "C" {
	//Cholesky decomposition of symmetric/hermitian positive definite matrix
	void spotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);
	void dpotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);
	void cpotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);
	void zpotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);


	//Inverse of triangular matrix
	void strtri_( char* UPLO, char* DIAG, int* N, void* A, int* LDA, int* INFO );
	void dtrtri_( char* UPLO, char* DIAG, int* N, void* A, int* LDA, int* INFO );
	void ctrtri_( char* UPLO, char* DIAG, int* N, void* A, int* LDA, int* INFO );
	void ztrtri_( char* UPLO, char* DIAG, int* N, void* A, int* LDA, int* INFO );
}

/**
 *   Perform Cholesky decomposition of matrix.
 *   hoNDArray should be symmetric/hermitian positive definite.
 *   Matrix will be replaced with lower triangular matrix.
 *   Calls LAPACK subroutine _POTRF
 *
 */
template <typename T> EXPORTLINALG int hoNDArray_choldc(hoNDArray<T>* A);

/**
 * Invert matrix assuming it is lower trinagular
 */
template <typename T> EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray<T>* A);

/**
 * Transpose matrix (if it is indeed a matrix (2 dimensional))
 */
template <typename T> EXPORTLINALG int hoNDArray_transpose(hoNDArray<T>* A_in, hoNDArray<T>* A_out);

#endif /* MATRIX_DECOMPOSITION_H_ */
