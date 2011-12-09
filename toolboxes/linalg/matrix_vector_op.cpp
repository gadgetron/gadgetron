/*
 * matrix_vector_op.cpp
 *
 *  Created on: Dec 9, 2011
 *      Author: Michael S. Hansen
 */

#include "matrix_vector_op.h"

	/*
	void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	                 const int K, const void *alpha, const void  *A,
	                 const int lda, const void  *B, const int ldb,
	                 const void *beta, void  *C, const int ldc)
	*/

void cblas_gemm_wrapper(std::complex<float>* A, std::complex<float>* B, std::complex<float>* C,
		int M, int N, int K, std::complex<float>* alpha, std::complex<float>* beta)
{
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, K, B, N, beta, C, N);
}

void cblas_gemm_wrapper(std::complex<double>* A, std::complex<double>* B, std::complex<double>* C,
		int M, int N, int K, std::complex<double>* alpha, std::complex<double>* beta)
{
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, K, B, N, beta, C, N);
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
	cblas_gemm_wrapper(A->get_data_ptr(), B->get_data_ptr(), C->get_data_ptr(), M, N, K, &alpha, &beta);

	return 0;
}

//Template instanciations
template int hoNDArray_gemm( hoNDArray< std::complex<float> >* A, hoNDArray< std::complex<float> >* B, std::complex<float> alpha,  hoNDArray< std::complex<float> >* C, std::complex<float> beta);
template int hoNDArray_gemm( hoNDArray< std::complex<double> >* A, hoNDArray< std::complex<double> >* B, std::complex<double> alpha,  hoNDArray< std::complex<double> >* C, std::complex<double> beta);
