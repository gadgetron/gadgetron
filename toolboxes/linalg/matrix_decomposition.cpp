/*
 * matrix_decomposition.cpp
 *
 *  Created on: Dec 10, 2011
 *      Author: Michael S. Hansen
 */

#include "matrix_decomposition.h"
#include <complex>

void potrf_wrapper(char* UPLO, int* N, float* A, int* LDA, int* info)
{
	spotrf_(UPLO, N, A, LDA, info);
}

void potrf_wrapper(char* UPLO, int* N, double* A, int* LDA, int* info)
{
	dpotrf_(UPLO, N, A, LDA, info);
}

void potrf_wrapper(char* UPLO, int* N, std::complex<float>* A, int* LDA, int* info)
{
	cpotrf_(UPLO, N, A, LDA, info);
}

void potrf_wrapper(char* UPLO, int* N, std::complex<double>* A, int* LDA, int* info)
{
	zpotrf_(UPLO, N, A, LDA, info);
}


template <typename T> int hoNDArray_choldc(hoNDArray<T>* A)
{
	char UPLO = 'L';
	int N = A->get_size(0);
	int LDA = N;
	int info = 0;

	potrf_wrapper(&UPLO, &N, A->get_data_ptr(), &LDA, &info);

	return info;
}


//Template instanciations
template int hoNDArray_choldc(hoNDArray< std::complex<float> >* A);
template int hoNDArray_choldc(hoNDArray< std::complex<double> >* A);
template int hoNDArray_choldc(hoNDArray< float >* A);
template int hoNDArray_choldc(hoNDArray< double >* A);




