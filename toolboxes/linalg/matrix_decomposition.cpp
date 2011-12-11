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
	const char* fname = "hoNDArray_choldc(hoNDArray<T>* A)";
	/*
	 *  We are specifying Upper Triangular,
	 *  but matrix comes in transposed (row-major) compared to
	 *  Fortran column-major order. As a result, we will get the lower
	 *  triangular matrix.
	 */
	char UPLO = 'U';
	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": This is not a matrix, only two dimensions allowed\n" << std::endl;
		return -1;
	}

	int N = A->get_size(0);
	if (N != A->get_size(1)) {
		std::cout << fname << ": Matrix is not symmetric.\n" << std::endl;
		return -1;
	}

	int LDA = N;
	int info = 0;

	potrf_wrapper(&UPLO, &N, A->get_data_ptr(), &LDA, &info);

	if (info != 0) {
		std::cout << fname << ": Error calling _potrf wrapper routine.\n" << std::endl;
		return -1;
	}

	return info;
}


//Template instanciations
template int hoNDArray_choldc(hoNDArray< std::complex<float> >* A);
template int hoNDArray_choldc(hoNDArray< std::complex<double> >* A);
template int hoNDArray_choldc(hoNDArray< float >* A);
template int hoNDArray_choldc(hoNDArray< double >* A);




