/*
 * matrix_vector_op.cpp
 *
 *  Created on: Dec 9, 2011
 *      Author: Michael S. Hansen
 */

#include "matrix_vector_op.h"
#include <cstdlib>

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
  
  //AXPY - Generalized matrix-matrix multiplication
  void saxpy_(int* N, void* A, void* X, int* INCX, void* Y, int* INCY);
  void daxpy_(int* N, void* A, void* X, int* INCX, void* Y, int* INCY);
  void caxpy_(int* N, void* A, void* X, int* INCX, void* Y, int* INCY);
  void zaxpy_(int* N, void* A, void* X, int* INCX, void* Y, int* INCY);
  
  //NORM2
  float snrm2_(int*N, void* X, int* INCX);
  float scnrm2_(int*N, void* X, int* INCX);
  double dnrm2_(int*N, void* X, int* INCX);
  double dznrm2_(int*N, void* X, int* INCX);

  //ASUM
  float sasum_(int* N, void* SX, int* INCX);
  float scasum_(int* N, void* SX, int* INCX);
  double dasum_(int* N, void* SX, int* INCX);
  double dzasum_(int* N, void* SX, int* INCX);

}

namespace Gadgetron
{
  double cblas_norm2_wrapper(int N, float* X, int INCX) {
    return static_cast<double>(snrm2_(&N,X,&INCX));
  }

  double cblas_norm2_wrapper(int N, double* X, int INCX) {
    return dnrm2_(&N,X,&INCX);
  }

  double  cblas_norm2_wrapper(int N, std::complex<float>* X, int INCX) {
    return static_cast<double>(scnrm2_(&N,X,&INCX));
  }

  double  cblas_norm2_wrapper(int N, std::complex<double>* X, int INCX) {
    return dznrm2_(&N,X,&INCX);
  }

  template <typename T> double hoNDArray_norm2(hoNDArray<T>* X)
  {
    return cblas_norm2_wrapper(X->get_number_of_elements(),X->get_data_ptr(),1); 
  }

  //Template instanciations
  template EXPORTLINALG double hoNDArray_norm2( hoNDArray<float>* X);
  template EXPORTLINALG double hoNDArray_norm2( hoNDArray<double>* X);
  template EXPORTLINALG double hoNDArray_norm2( hoNDArray< std::complex<float> >* X);
  template EXPORTLINALG double hoNDArray_norm2( hoNDArray< std::complex<double> >* X);

  
  double cblas_asum_wrapper(int N, float* X, int INCX) {
    return static_cast<double>(sasum_(&N,X,&INCX));
  }

  double cblas_asum_wrapper(int N, double* X, int INCX) {
    return dasum_(&N,X,&INCX);
  }

  double  cblas_asum_wrapper(int N, std::complex<float>* X, int INCX) {
    return static_cast<double>(scasum_(&N,X,&INCX));
  }

  double  cblas_asum_wrapper(int N, std::complex<double>* X, int INCX) {
    return dzasum_(&N,X,&INCX);
  }

 template <typename T> double hoNDArray_asum(hoNDArray<T>* X)
  {
    return cblas_asum_wrapper(X->get_number_of_elements(),X->get_data_ptr(),1); 
  }

  //Template instanciations
  template EXPORTLINALG double hoNDArray_asum( hoNDArray<float>* X);
  template EXPORTLINALG double hoNDArray_asum( hoNDArray<double>* X);
  template EXPORTLINALG double hoNDArray_asum( hoNDArray< std::complex<float> >* X);
  template EXPORTLINALG double hoNDArray_asum( hoNDArray< std::complex<double> >* X);

  void cblas_axpy_wrapper(int N, float* A, float* X, int INCX, float* Y, int INCY)
  {
    saxpy_(&N, A, X, &INCX, Y, &INCY);
  }

  void cblas_axpy_wrapper(int N, double* A, double* X, int INCX, double* Y, int INCY)
  {
    daxpy_(&N, A, X, &INCX, Y, &INCY);
  }

  void cblas_axpy_wrapper(int N, std::complex<float>* A, std::complex<float>* X, int INCX, std::complex<float>* Y, int INCY)
  {
    caxpy_(&N, A, X, &INCX, Y, &INCY);
  }

  void cblas_axpy_wrapper(int N, std::complex<double>* A, std::complex<double>* X, int INCX, std::complex<double>* Y, int INCY)
  {
    zaxpy_(&N, A, X, &INCX, Y, &INCY);
  }

  template <typename T> void hoNDArray_axpy( T* A, hoNDArray<T>* X, hoNDArray<T>* Y)
  {
    cblas_axpy_wrapper((int)X->get_number_of_elements(),A,X->get_data_ptr(),1,Y->get_data_ptr(),1);
  }

  //Template instanciations
  template EXPORTLINALG void hoNDArray_axpy( float* A, hoNDArray< float >* X, hoNDArray< float>* Y);
  template EXPORTLINALG void hoNDArray_axpy( double* A, hoNDArray< double >* X, hoNDArray< double>* Y);

  template EXPORTLINALG void hoNDArray_axpy( std::complex<float>* A, hoNDArray< std::complex<float> >* X, 
					     hoNDArray< std::complex<float> >* Y);

  template EXPORTLINALG void hoNDArray_axpy( std::complex<double>* A, hoNDArray< std::complex<double> >* X, 
					     hoNDArray< std::complex<double> >* Y);

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
  
  template <typename T> void hoNDArray_gemm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha,  hoNDArray<T>* C, T beta)
  {

	//Let's first check the dimensions A
	if (A->get_number_of_dimensions() != 2) {
		throw std::runtime_error("Invalid number of dimensions in matrix A: ");
	}

	//Let's first check the dimensions B
	if (B->get_number_of_dimensions() != 2) {
		throw std::runtime_error("Invalid number of dimensions in matrix B: ");
	}

	//Let's first check the dimensions C
	if (C->get_number_of_dimensions() != 2) {
		throw std::runtime_error("Invalid number of dimensions in matrix C: ");
	}

	//Do the dimensions match?
	int M = A->get_size(1); //Number of rows of A
	int N = B->get_size(0); //Number of columns of B
	int K = A->get_size(0); //Number of columns of A

	if (K != static_cast<int>(B->get_size(1))) {
		std::stringstream ss;
		ss <<"Number of columns of A (" << K << ") does not match rows of B(" << B->get_size(1) << ")" << std::endl;
		throw std::runtime_error(ss.str());
	}


	//Is the output matric the right size?
	if ((C->get_size(0) != N) || (C->get_size(1) != M) ) {
		std::stringstream ss;
		ss << "Size of output matrix C (" << C->get_size(0) << " (cols), " <<
		C->get_size(1) << " (rows))" << " does not match the expected" <<
		N << "(cols), " << M << "(rows)" << std::endl;
		throw std::runtime_error(ss.str());


	}

	//Now call appropriate CBLAS routine
	char TRANSA = 'N';
	char TRANSB = 'N';
	cblas_gemm_wrapper(TRANSA, TRANSB, A->get_data_ptr(), B->get_data_ptr(), C->get_data_ptr(), M, N, K, &alpha, &beta);

}

//Template instanciations
template EXPORTLINALG void hoNDArray_gemm( hoNDArray< float>* A, hoNDArray< float >* B, float alpha,  hoNDArray< float >* C, float beta);
template EXPORTLINALG void hoNDArray_gemm( hoNDArray< double >* A, hoNDArray< double >* B, double alpha,  hoNDArray< double >* C, double beta);
template EXPORTLINALG void hoNDArray_gemm( hoNDArray< std::complex<float> >* A, hoNDArray< std::complex<float> >* B, std::complex<float> alpha,  hoNDArray< std::complex<float> >* C, std::complex<float> beta);
template EXPORTLINALG void hoNDArray_gemm( hoNDArray< std::complex<double> >* A, hoNDArray< std::complex<double> >* B, std::complex<double> alpha,  hoNDArray< std::complex<double> >* C, std::complex<double> beta);

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

template <typename T> void hoNDArray_trmm( hoNDArray<T>* A, hoNDArray<T>* B, T alpha)
{
	//Let's first check the dimensions A
	if (A->get_number_of_dimensions() != 2) {
		throw std::runtime_error("Invalid number of dimensions in matrix A: ");
	}

	//Let's first check the dimensions B
	if (B->get_number_of_dimensions() != 2) {
		throw std::runtime_error("Invalid number of dimensions in matrix B: ");
	}

	//Do the dimensions match?
	int M = A->get_size(1); //Number of rows of A
	int N = B->get_size(0); //Number of columns of B

	trmm_wrapper(&M, &N, &alpha, A->get_data_ptr(), B->get_data_ptr());

	}

template EXPORTLINALG void hoNDArray_trmm( hoNDArray<float>* A, hoNDArray<float>* B, float alpha);
template EXPORTLINALG void hoNDArray_trmm( hoNDArray<double>* A, hoNDArray<double>* B, double alpha);
template EXPORTLINALG void hoNDArray_trmm( hoNDArray< std::complex<float> >* A, hoNDArray< std::complex<float> >* B, std::complex<float> alpha);
template EXPORTLINALG void hoNDArray_trmm( hoNDArray< std::complex<double> >* A, hoNDArray< std::complex<double> >* B, std::complex<double> alpha);

} //Namespace Gadgetron
