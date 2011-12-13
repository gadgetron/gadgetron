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

	/* Temp code to zero upper triangular */
	/*
	T* d = A->get_data_ptr();
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = i+1; j < N; j++) {
			d[i*N+j] = 0;
		}
	}
	*/

	return info;
}


//Template instanciations
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< std::complex<float> >* A);
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< std::complex<double> >* A);
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< float >* A);
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< double >* A);


void trtri_wrapper(char* UPLO, char* DIAG, int* N, float* A, int* LDA, int* info)
{
	strtri_(UPLO, DIAG, N, A, LDA, info);
}

void trtri_wrapper(char* UPLO, char* DIAG, int* N, double* A, int* LDA, int* info)
{
	dtrtri_(UPLO, DIAG, N, A, LDA, info);
}

void trtri_wrapper(char* UPLO, char* DIAG, int* N, std::complex<float>* A, int* LDA, int* info)
{
	ctrtri_(UPLO, DIAG, N, A, LDA, info);
}

void trtri_wrapper(char* UPLO, char* DIAG, int* N, std::complex<double>* A, int* LDA, int* info)
{
	ztrtri_(UPLO, DIAG, N, A, LDA, info);
}

template <typename T> int hoNDArray_inv_lower_triangular(hoNDArray<T>* A)
{
	const char* fname = "hoNDArray_inv_lower_triangular(hoNDArray<T>* A)";

	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": Error array is not 2 dimensional.\n" << std::endl;
		return -1;
	}

	int N = A->get_size(0);

	if (N != A->get_size(1)) {
		std::cout << fname << ": Error array is not 2 dimensional.\n" << std::endl;
		return -1;
	}

	int LDA = N;
	char UPLO = 'U'; //We are passing upper, but matrix is really lower. This is do deal with row and column major order differences
	char DIAG = 'N';
	int info;

	trtri_wrapper(&UPLO, &DIAG, &N, A->get_data_ptr(), &LDA, &info);

	if (info != 0) {
		std::cout << fname << ": Error inverting triangular matrix.\n" << std::endl;
		return -1;
	}

	return 0;
}

template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray<float>* A);
template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray<double>* A);
template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray< std::complex<float> >* A);
template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray< std::complex<double> >* A);


template<typename T> int hoNDArray_svd(hoNDArray<T> *A, hoNDArray<T> *U, hoNDArray<T> *S, hoNDArray<T> *VT)
{
	const char* fname = "hoNDArray_svd(hoNDArray<T> *A, hoNDArray<T> *U, hoNDArray<T> *S, hoNDArray<T> *VT): ";

	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": Error array is not 2 dimensional.\n" << std::endl;
		return -1;
	}

	int M = A->get_size(1);
	int N = A->get_size(0);

	std::vector<unsigned int> permute_dims(2);
	permute_dims[0] = 1;permute_dims[1] = 0;

	std::vector<unsigned int> perm_dims(2);
	perm_dims[0] = M;
	perm_dims[1] = N;

	hoNDArray< T > A_perm;
	if (!A_perm.create(&perm_dims)) {
		std::cout << fname << ": Unable to create storage for transposed array\n" << std::endl;
		return -1;
	}

	A->permute(&permute_dims, &A_perm);
	T* A_ptr = A_perm.get_data_ptr();

	char JOBU, JOBVT;
	T* U_ptr = 0;
	T* VT_ptr = 0;

	int U_M = 0;
	int U_N = 0;

	if (U) {
		U_M = U->get_size(1);
		U_N = U->get_size(0);

		//if ((U_M == U_N) && )

	} else {
		JOBU = 'N';
	}

	int VT_M = 0;
	int VT_N = 0;
	if (VT) {

	} else {
		JOBVT = 'N';
	}

	return 0;
}


template int hoNDArray_svd(hoNDArray< std::complex<float> > *A, hoNDArray< std::complex<float> > *U, hoNDArray< std::complex<float> > *S, hoNDArray< std::complex<float> > *VT);

