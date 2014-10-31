
#pragma once

#include "cpucore_math_export.h"
#include "hoMatrix.h"

#ifdef USE_ARMADILLO
    #include "hoArmadillo.h"
#endif // USE_ARMADILLO

#ifndef lapack_int
    #define lapack_int int
#endif // lapack_int

/// ----------------------------------------------------------------------
/// the fortran interface of lapack and blas functions are called
/// ----------------------------------------------------------------------

namespace Gadgetron
{

template<typename T> EXPORTCPUCOREMATH
bool EigenAnalysis_syev_heev2(hoMatrix<T>& A, hoMatrix<T>& eigenValue);

template<typename T> EXPORTCPUCOREMATH
bool SolveLinearSystem_Tikhonov(hoMatrix<T>& A, hoMatrix<T>& b, hoMatrix<T>& x, double lamda);

// ----------------------------------------------------------------------------

// following matrix computation calls lapacke functions

/// C = A*B for complex float
EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm_CXFL(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B);

template<typename T> EXPORTCPUCOREMATH
bool GeneralMatrixProduct_gemm(hoNDArray<T>& C, 
                            const hoNDArray<T>& A, bool transA, 
                            const hoNDArray<T>& B, bool transB);

//template<typename T> EXPORTCPUCOREMATH 
//bool GeneralMatrixProduct_gemm(hoMatrix<T>& C, 
//                            const hoMatrix<T>& A, bool transA, 
//                            const hoMatrix<T>& B, bool transB);

/// Performs a symmetric rank-k update (no conjugated).
template<typename T> EXPORTCPUCOREMATH 
bool syrk(hoNDArray<T>& C, const hoNDArray<T>& A, char uplo, bool isATA);

/// Performs a Hermitian rank-k update.
template<typename T> EXPORTCPUCOREMATH 
bool herk(hoNDArray<T>& C, const hoNDArray<T>& A, char uplo, bool isAHA);

template<typename T> EXPORTCPUCOREMATH 
bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo);

template<typename T> EXPORTCPUCOREMATH 
bool EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue);

template<typename T> EXPORTCPUCOREMATH 
bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A);

template<typename T> EXPORTCPUCOREMATH 
bool TriangularInverse_trtri(hoMatrix<T>& A, char uplo);

template<typename T> EXPORTCPUCOREMATH 
bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b);

/// Computes the LU factorization of a general m-by-n matrix
/// this function is called by general matrix inversion
template<typename T> EXPORTCPUCOREMATH 
bool LUFactorizationGeneralMatrix_getrf(hoMatrix<T>& A, hoNDArray<lapack_int>& ipiv);

/// Computes the inverse of an LU-factored general matrix
template<typename T> EXPORTCPUCOREMATH 
bool InverseGeneralMatrix_getri(hoMatrix<T>& A);

}
