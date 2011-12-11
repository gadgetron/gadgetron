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

//Declaration of lapack routine
extern "C" {
	//Cholesky decomposition of symmetric/hermitian positive definite matrix
	void spotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);
	void dpotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);
	void cpotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);
	void zpotrf_(char* UPLO, int* N, void* A, int* LDA, int* info);
}

/**
 *   Perform Cholesky decomposition of matrix.
 *   hoNDArray should be symmetric/hermitian positive definite.
 *   Matrix will be replaced with lower triangular matrix.
 *   Calls LAPACK subroutine _POTRF
 *
 */
template <typename T> EXPORTLINALG int hoNDArray_choldc(hoNDArray<T>* A);

#endif /* MATRIX_DECOMPOSITION_H_ */
