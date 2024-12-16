
#pragma once

#include "hoNDArray.h"

#include "cpp_lapack.h"
//#include "hoArmadillo.h"



/// ----------------------------------------------------------------------
/// the fortran interface of lapack and blas functions are called
/// ----------------------------------------------------------------------

namespace Gadgetron
{

// following matrix computation calls lapacke functions

/// C = A*B
template<class T>
void gemm(hoNDArray<T>& C, const hoNDArray<T>& A, const hoNDArray<T>& B);
/// if transA==true, C = A'*B
/// if transB==true, C=A*B'
/// if both are true, C=A'*B'
template<typename T>
void gemm(hoNDArray<T>& C,
        const hoNDArray<T>& A, bool transA,
        const hoNDArray<T>& B, bool transB);

/// perform a symmetric rank-k update (no conjugated).
template<typename T>
void syrk(hoNDArray<T>& C, const hoNDArray<T>& A, char uplo, bool isATA);

/// perform a Hermitian rank-k update.
template<typename T>
void herk(hoNDArray<T>& C, const hoNDArray<T>& A, char uplo, bool isAHA);

/// compute the Cholesky factorization of a real symmetric positive definite matrix A
template<typename T>
void potrf(hoNDArray<T>& A, char uplo);

/// compute all eigenvalues and eigenvectors of a Hermitian matrix A
template<typename T>
std::enable_if_t<std::is_floating_point_v<T>> heev(hoNDArray<T>& A, hoNDArray<T>& eigenValue);

template<typename T>
void heev(hoNDArray< std::complex<T> >& A, hoNDArray<T>& eigenValue);

/// compute inverse of a symmetric (Hermitian) positive-definite matrix A
template<typename T>
void potri(hoNDArray<T>& A);

/// compute the inverse of a triangular matrix A
template<typename T>
void trtri(hoNDArray<T>& A, char uplo);

/// solve Ax=b, a symmetric or Hermitian positive-definite matrix A and multiple right-hand sides b
/// b is replaced with x
template<typename T>
void posv(hoNDArray<T>& A, hoNDArray<T>& b);

/// solve Ax=b, a square symmetric / hermitian matrix A and multiple right-hand sides b
/// for float and double, A is a symmetric matrix
/// for complex type, A is a hermitian matrix
/// b is replaced with x
template<typename T>
void hesv(hoNDArray<T>& A, hoNDArray<T>& b);

/// solve Ax=b, a square matrix A and multiple right-hand sides b
/// b is replaced with x
template<typename T>
void gesv(hoNDArray<T>& A, hoNDArray<T>& b);

/// solve Ax=b with Tikhonov regularization
template<typename T>
void SolveLinearSystem_Tikhonov(hoNDArray<T>& A, hoNDArray<T>& b, hoNDArray<T>& x, double lamda);

/// Computes the LU factorization of a general m-by-n matrix
/// this function is called by general matrix inversion
template<typename T>
void getrf(hoNDArray<T>& A, hoNDArray<Lapack::Int>& ipiv);

/// Computes the inverse of a matrix using LU factorization
template<typename T>
void invert(hoNDArray<T>& A);

/**
* @brief linear fitting, y = a*x + b
  compute linear fit for y to x
*/
template <typename T>  void linFit(const hoNDArray<T>& x, const hoNDArray<T>& y, T& a, T& b);

/**
* @brief 2nd-order poly fitting, y = a*x^2 + b*x + c
  compute linear fit for y to x
*/
template <typename T>  void polyFit2(const hoNDArray<T>& x, const hoNDArray<T>& y, T& a, T& b, T& c);

}
