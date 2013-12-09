#pragma once

#include "ho2DArray.h"
#include "complext.h"

#ifdef USE_MKL
    #include "mkl.h"
#endif // USE_MKL

#ifdef GT_Complex8
    #undef GT_Complex8
#endif // GT_Complex8
typedef std::complex<float> GT_Complex8;

#ifdef GT_Complex16
    #undef GT_Complex16
#endif // GT_Complex16
typedef std::complex<double> GT_Complex16; 

namespace Gadgetron{

// the hoMatrix stores every row as the first dimension
template <class T> class  hoMatrix : public ho2DArray<T>
{
public:

    typedef hoMatrix<T> Self;
    typedef ho2DArray<T> BaseClass;

    hoMatrix();
    hoMatrix(unsigned long long rows, unsigned long long cols);
    hoMatrix(unsigned long long rows, unsigned long long cols, T* data, bool delete_data_on_destruct = false);

    virtual ~hoMatrix();

    hoMatrix(const hoMatrix<T>& a);
    hoMatrix<T>& operator=(const hoMatrix& rhs);

    virtual bool createMatrix(unsigned long long rows, unsigned long long cols);
    virtual bool createMatrix(unsigned long long rows, unsigned long long cols, T* data, bool delete_data_on_destruct = false);

    T& operator()(long long r , long long c);
    const T& operator()(long long r , long long c) const;

    unsigned long long rows() const;
    unsigned long long cols() const;

    // assign the upper/lower triangle matrix as a fixed value
    bool upperTri(const T& v);
    bool lowerTri(const T& v);

    // sum along row or col
    bool sumOverRow(hoNDArray<T>& res) const;
    bool sumOverCol(hoNDArray<T>& res) const;

    // get the sub matrix
    bool subMatrix(Self& res, unsigned long long startR, unsigned long long endR, unsigned long long startC, unsigned long long endC) const;

    bool operator == (const Self& m) const;
    bool operator != (const Self& m) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;
    using BaseClass::accesser_;
};

// copy lower half of a matrix to its upper half
//template<typename T>  bool copyL2U(hoMatrix<T>& A);
//// conj: true, perform conjuate while copying
//template<typename T>  bool copyL2U(hoMatrix<T>& A, bool conj);
//
//// copy upper half of a matrix to its lower half
//template<typename T>  bool copyU2L(hoMatrix<T>& A);
//template<typename T>  bool copyU2L(hoMatrix<T>& A, bool conj);
//
//// matrix transpose/conjugate transpose
//template<typename T>  bool trans(const hoMatrix<T>& A, hoMatrix<T>& AT);
//template<typename T>  bool conjugatetrans(const hoMatrix<T>& A, hoMatrix<T>& AH);
//template<>  bool conjugatetrans(const hoMatrix<float>& A, hoMatrix<float>& AH);
//template<>  bool conjugatetrans(const hoMatrix<double>& A, hoMatrix<double>& AH);
//
//// following matrix computation calls MKL functions
//#ifdef USE_MKL
//
//    // All description below are from MKL manual!
//
//    // BLAS gemm : matrix-matrix product with general matrices
//    // The gemm routines compute a scalar-matrix-matrix product and add the result to a scalar-matrix product, with general matrices. The operation is defined as
//    // C := alpha*op(A)*op(B) + beta*C,
//    // where:
//    // op(X) is one of op(X) = X, or op(X) = XT, or op(X) = XH,
//    // alpha and beta are scalars,
//    // A, B and C are matrices:
//    // op(A) is an m-by-k matrix,
//    // op(B) is a k-by-n matrix,
//    // C is an m-by-n matrix.
//    template<typename T>  bool GeneralMatrixProduct_gemm(hoMatrix<T>& C, const hoMatrix<T>& A, bool transA, const hoMatrix<T>& B, bool transB);
//
//    // LAPACK potrf: Cholesky factorization of a symmetric (Hermitian) positive-definite matrix
//    // description from MKL manual:
//    // The routine forms the Cholesky factorization of a symmetric positive-definite or, for complex data, Hermitian positive-definite matrix A:
//    // A = UT*U for real data, A = UH*U for complex data, if uplo='U'
//    // A = L*LT for real data, A = L*LH for complex data, if uplo='L'
//    // where L is a lower triangular matrix and U is upper triangular.
//    // on return, the upper or lower triangular part of A is overwritten by the Cholesky factor U or L, as specified by uplo.
//    template<typename T>  bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo);
//
//    // LAPACK syev/heev
//    // Computes all eigenvalues and, optionally, eigenvectors of a real symmetric/complex hermitian matrix.
//    // on return, A stores the eigen vectors as column vectors, D stores the eigen values
//    template<typename T>  bool EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue);
//    template<typename T>  bool EigenAnalysis_syev_heev2(hoMatrix<T>& A, hoMatrix<T>& eigenValue);
//
//    // LAPACK potri
//    // Computes the inverse of a symmetric (Hermitian) positive-definite matrix.
//    // on return, A holds the inverse matrix
//    template<typename T>  bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A);
//
//    // LAPACK trtri
//    // Computes the inverse of a triangular matrix.
//    // on return, A holds the inverse matrix
//    template<typename T>  bool TriangularInverse_trtri(hoMatrix<T>& A, char uplo);
//
//    // LAPACK posv
//    // Computes the solution to the system of linear equations with a symmetric or Hermitian positive-definite matrix A and multiple right-hand sides.
//    // Ax=b
//    // on return, A is replaced by its Cholesky factorization and b holds the solution
//    template<typename T>  bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b);
//
//    // Lineare equation solve with regularization
//    // Ax = b
//    // what is solved is (A'A+lamda*I)x = A'b
//    template<typename T>  bool SolveLinearSystem_Tikhonov(hoMatrix<T>& A, hoMatrix<T>& b, hoMatrix<T>& x, double lamda);
//
//#endif // USE_MKL

}

#include <hoMatrix.cpp>
