/*
 * matrix_vector_op.cpp
 *
 *  Created on: Dec 9, 2011
 *      Author: Michael S. Hansen
 */

#include "matrix_vector_op.h"

void cblas_gemm_wrapper(char TRANSA, char TRANSB, float* A, float* B, float* C,
		int M, int N, int K, float* alpha, float* beta)
{
	//We have to flip the arguments here to make it fit with FORTRAN column-major order so that we don't have to transpose
	sgemm_(&TRANSB, &TRANSA,&N, &M, &K, alpha, B, &N, A, &K, beta, C, &N);
}

void cblas_gemm_wrapper(char TRANSA, char TRANSB, double* A, double* B, double* C,
		int M, int N, int K, double* alpha, double* beta)
{
	//We have to flip the arguments here to make it fit with FORTRAN column-major order so that we don't have to transpose
	dgemm_(&TRANSB, &TRANSA,&N, &M, &K, alpha, B, &N, A, &K, beta, C, &N);
}

void cblas_gemm_wrapper(char TRANSA, char TRANSB, std::complex<float>* A, std::complex<float>* B, std::complex<float>* C,
		int M, int N, int K, std::complex<float>* alpha, std::complex<float>* beta)
{
	//We have to flip the arguments here to make it fit with FORTRAN column-major order so that we don't have to transpose
	cgemm_(&TRANSB, &TRANSA,&N, &M, &K, alpha, B, &N, A, &K, beta, C, &N);
}

void cblas_gemm_wrapper(char TRANSA, char TRANSB, std::complex<double>* A, std::complex<double>* B, std::complex<double>* C,
		int M, int N, int K, std::complex<double>* alpha, std::complex<double>* beta)
{
	//We have to flip the arguments here to make it fit with FORTRAN column-major order so that we don't have to transpose
	zgemm_(&TRANSB, &TRANSA,&N, &M, &K, alpha, B, &N, A, &K, beta, C, &N);
}

template <typename T> int hoNDArray_gemm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha,  hoNDArray<T>* C, T beta)
{
	const char* fname = "hoNDArray_gemm";

	//Let's first check the dimensions A
	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": " << "Invalid number of dimensions in matrix A: " << A->get_number_of_dimensions() << std::endl;
		return -1;
	}

	//Let's first check the dimensions B
	if (B->get_number_of_dimensions() != 2) {
		std::cout << fname << ": " << "Invalid number of dimensions in matrix B: " << B->get_number_of_dimensions() << std::endl;
		return -1;
	}

	//Let's first check the dimensions C
	if (C->get_number_of_dimensions() != 2) {
		std::cout << fname << ": " << "Invalid number of dimensions in matrix C: " << C->get_number_of_dimensions() << std::endl;
		return -1;
	}

	//Do the dimensions match?
	int M = A->get_size(1); //Number of rows of A
	int N = B->get_size(0); //Number of columns of B
	int K = A->get_size(0); //Number of columns of A

	if (K != static_cast<int>(B->get_size(1))) {
		std::cout << fname << ": " << "Number of columns of A (" << K << ") does not match rows of B("
				<< B->get_size(1) << ")" << std::endl;
		return -1;
	}


	//Is the output matric the right size?
	if ((C->get_size(0) != N) || (C->get_size(1) != M) ) {
		std::cout << fname << ": " << "Size of output matrix C (" << C->get_size(0) << " (cols), "
				<< C->get_size(1) << " (rows))" << " does not match the expected" <<
				N << "(cols), " << M << "(rows)" << std::endl;
		return -1;

	}

	//Now call appropriate CBLAS routine
	char TRANSA = 'N';
	char TRANSB = 'N';
	cblas_gemm_wrapper(TRANSA, TRANSB, A->get_data_ptr(), B->get_data_ptr(), C->get_data_ptr(), M, N, K, &alpha, &beta);

	return 0;
}

//Template instanciations
template EXPORTLINALG int hoNDArray_gemm( hoNDArray< float>* A, hoNDArray< float >* B, float alpha,  hoNDArray< float >* C, float beta);
template EXPORTLINALG int hoNDArray_gemm( hoNDArray< double >* A, hoNDArray< double >* B, double alpha,  hoNDArray< double >* C, double beta);
template EXPORTLINALG int hoNDArray_gemm( hoNDArray< std::complex<float> >* A, hoNDArray< std::complex<float> >* B, std::complex<float> alpha,  hoNDArray< std::complex<float> >* C, std::complex<float> beta);
template EXPORTLINALG int hoNDArray_gemm( hoNDArray< std::complex<double> >* A, hoNDArray< std::complex<double> >* B, std::complex<double> alpha,  hoNDArray< std::complex<double> >* C, std::complex<double> beta);

void trmm_wrapper(int* M,int* N, float* ALPHA,float* A, void* B) {
	char SIDE = 'R'; char UPLO = 'U'; char TRANSA = 'N'; char DIAG = 'N';
	strmm_(&SIDE, &UPLO, &TRANSA, &DIAG, N, M, ALPHA, B, N, A, M);
}

void trmm_wrapper(int* M,int* N, double* ALPHA, double* A, void* B) {
	char SIDE = 'R'; char UPLO = 'U'; char TRANSA = 'N'; char DIAG = 'N';
	dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, N, M, ALPHA, B, N, A, M);
}

void trmm_wrapper(int* M,int* N, std::complex<float>* ALPHA,std::complex<float>* A, void* B) {
	char SIDE = 'R'; char UPLO = 'U'; char TRANSA = 'N'; char DIAG = 'N';
	ctrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, N, M, ALPHA, A, M, B, N);
}

void trmm_wrapper(int* M,int* N, std::complex<double>* ALPHA,std::complex<double>* A, void* B) {
	char SIDE = 'R'; char UPLO = 'U'; char TRANSA = 'N'; char DIAG = 'N';
	ztrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, N, M, ALPHA, B, N, A, M);
}

template <typename T> int hoNDArray_trmm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha)
{
	const char* fname = "hoNDArray_trmm";

	//Let's first check the dimensions A
	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": " << "Invalid number of dimensions in matrix A: " << A->get_number_of_dimensions() << std::endl;
		return -1;
	}

	//Let's first check the dimensions B
	if (B->get_number_of_dimensions() != 2) {
		std::cout << fname << ": " << "Invalid number of dimensions in matrix B: " << B->get_number_of_dimensions() << std::endl;
		return -1;
	}

	//Do the dimensions match?
	int M = A->get_size(1); //Number of rows of A
	int N = B->get_size(0); //Number of columns of B

	trmm_wrapper(&M, &N, &alpha, A->get_data_ptr(), B->get_data_ptr());

	return 0;
}

template EXPORTLINALG int hoNDArray_trmm( hoNDArray<float>* A, hoNDArray<float>* B, float alpha);
template EXPORTLINALG int hoNDArray_trmm( hoNDArray<double>* A, hoNDArray<double>* B, double alpha);
template EXPORTLINALG int hoNDArray_trmm( hoNDArray< std::complex<float> >* A, hoNDArray< std::complex<float> >* B, std::complex<float> alpha);
template EXPORTLINALG int hoNDArray_trmm( hoNDArray< std::complex<double> >* A, hoNDArray< std::complex<double> >* B, std::complex<double> alpha);
