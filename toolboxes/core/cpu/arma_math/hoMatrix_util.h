
#pragma once

#include "cpucore_math_export.h"
#include "hoMatrix.h"

#ifdef USE_MKL
    #include "mkl.h"
#endif // USE_MKL

#ifdef USE_ARMADILLO
#include "hoArmadillo.h"
#endif // USE_ARMADILLO

#ifdef GT_Complex8
#undef GT_Complex8
#endif // GT_Complex8
typedef std::complex<float> GT_Complex8;

#ifdef GT_Complex16
#undef GT_Complex16
#endif // GT_Complex16
typedef std::complex<double> GT_Complex16;

namespace Gadgetron
{

#if defined(USE_MKL) || defined(USE_ARMADILLO)

template<typename T> EXPORTCPUCOREMATH
bool EigenAnalysis_syev_heev2(hoMatrix<T>& A, hoMatrix<T>& eigenValue);

template<typename T> EXPORTCPUCOREMATH
bool SolveLinearSystem_Tikhonov(hoMatrix<T>& A, hoMatrix<T>& b, hoMatrix<T>& x, double lamda);

#endif // defined(USE_MKL) || defined(USE_ARMADILLO)

// following matrix computation calls MKL functions
#ifdef USE_MKL

/// C = A*B for complex float
EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm_CXFL(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B);

template<typename T> EXPORTCPUCOREMATH
bool GeneralMatrixProduct_gemm(hoNDArray<T>& C, 
                            const hoNDArray<T>& A, bool transA, 
                            const hoNDArray<T>& B, bool transB);

template<typename T> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoMatrix<T>& C, 
                            const hoMatrix<T>& A, bool transA, 
                            const hoMatrix<T>& B, bool transB);

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

#else

    // matrix computation calls armadillo
    #ifdef USE_ARMADILLO

    template<typename T> EXPORTCPUCOREMATH 
    bool GeneralMatrixProduct_gemm(hoNDArray<T>& C, 
                                const hoNDArray<T>& A, bool transA, 
                                const hoNDArray<T>& B, bool transB);

    template<typename T> EXPORTCPUCOREMATH 
    bool GeneralMatrixProduct_gemm(hoMatrix<T>& C, 
                                const hoMatrix<T>& A, bool transA, 
                                const hoMatrix<T>& B, bool transB);

    template<typename T> EXPORTCPUCOREMATH 
    bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A);

    template<typename T> EXPORTCPUCOREMATH 
    bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo);

    template<typename T> EXPORTCPUCOREMATH 
    bool EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue);

    template<typename T> EXPORTCPUCOREMATH 
    bool TriangularInverse_trtri(hoMatrix<T>& A, char uplo);

    template<typename T> EXPORTCPUCOREMATH 
    bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b);

    template<typename T> EXPORTCPUCOREMATH 
    bool InverseGeneralMatrix_getri(hoMatrix<T>& A);

    #endif // USE_ARMADILLO

#endif // USE_MKL

}
