
#include "hoMatrix_util.h"

#ifdef USE_MKL
    #include "mkl.h"
#endif // USE_MKL

#ifdef USE_MKL

    #define sgemm_ sgemm
    #define dgemm_ dgemm
    #define cgemm_ cgemm
    #define zgemm_ zgemm

#else
    #ifdef USE_LAPACK

        #ifndef lapack_complex_float
            #define lapack_complex_float GT_Complex8
        #endif // lapack_complex_float

        #ifndef lapack_complex_double
            #define lapack_complex_double GT_Complex16
        #endif // #ifndef lapack_complex_double

        extern "C" void sgemm_(const char *transa, const char *transb, const lapack_int *m, const lapack_int *n, const lapack_int *k,
                    const float *alpha, const float *a, const lapack_int *lda, const float *b, const lapack_int *ldb,
                    const float *beta, float *c, const lapack_int *ldc);

        extern "C" void dgemm_(const char *transa, const char *transb, const lapack_int *m, const lapack_int *n, const lapack_int *k,
                    const double *alpha, const double *a, const lapack_int *lda, const double *b, const lapack_int *ldb,
                    const double *beta, double *c, const lapack_int *ldc);

        extern "C" void cgemm_(const char *transa, const char *transb, const lapack_int *m, const lapack_int *n, const lapack_int *k,
                            const lapack_complex_float *alpha, const lapack_complex_float *a, const lapack_int *lda,
                            const lapack_complex_float *b, const lapack_int *ldb, const lapack_complex_float *beta,
                            lapack_complex_float *c, const lapack_int *ldc);

        extern "C" void zgemm_(const char *transa, const char *transb, const lapack_int *m, const lapack_int *n, const lapack_int *k,
                    const lapack_complex_double *alpha, const lapack_complex_double *a, const lapack_int *lda,
                    const lapack_complex_double *b, const lapack_int *ldb, const lapack_complex_double *beta,
                    lapack_complex_double *c, const lapack_int *ldc);

        extern "C" void spotrf_( const char* uplo, const lapack_int* n, float* a, const lapack_int* lda, lapack_int* info );
        extern "C" void dpotrf_( const char* uplo, const lapack_int* n, double* a, const lapack_int* lda, lapack_int* info );
        extern "C" void cpotrf_( const char* uplo, const lapack_int* n, lapack_complex_float* a, const lapack_int* lda, lapack_int* info );
        extern "C" void zpotrf_( const char* uplo, const lapack_int* n, lapack_complex_double* a, const lapack_int* lda, lapack_int* info );

        extern "C" void ssyev_( const char* jobz, const char* uplo, const lapack_int* n, float* a,
                const lapack_int* lda, float* w, float* work, const lapack_int* lwork,
                lapack_int* info );

        extern "C" void dsyev_( const char* jobz, const char* uplo, const lapack_int* n, double* a,
                const lapack_int* lda, double* w, double* work, const lapack_int* lwork,
                lapack_int* info );

        extern "C" void cheev_( const char* jobz, const char* uplo, const lapack_int* n,
                lapack_complex_float* a, const lapack_int* lda, float* w, lapack_complex_float* work,
                const lapack_int* lwork, float* rwork, lapack_int* info );

        extern "C" void zheev_( const char* jobz, const char* uplo, const lapack_int* n,
                lapack_complex_double* a, const lapack_int* lda, double* w,
                lapack_complex_double* work, const lapack_int* lwork, double* rwork,
                lapack_int* info );

        extern "C" void spotrf_( const char* uplo, const lapack_int* n, float* a, const lapack_int* lda,
                lapack_int* info );

        extern "C" void spotri_( const char* uplo, const lapack_int* n, float* a, const lapack_int* lda,
                lapack_int* info );

        extern "C" void dpotrf_( const char* uplo, const lapack_int* n, double* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void dpotri_( const char* uplo, const lapack_int* n, double* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void cpotrf_( const char* uplo, const lapack_int* n, lapack_complex_float* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void cpotri_( const char* uplo, const lapack_int* n, lapack_complex_float* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void zpotrf_( const char* uplo, const lapack_int* n, lapack_complex_double* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void zpotri_( const char* uplo, const lapack_int* n, lapack_complex_double* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void strtri_( const char* uplo, const char* diag, const lapack_int* n, float* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void dtrtri_( const char* uplo, const char* diag, const lapack_int* n, double* a,
                const lapack_int* lda, lapack_int* info );

        extern "C" void ctrtri_( const char* uplo, const char* diag, const lapack_int* n,
                lapack_complex_float* a, const lapack_int* lda, lapack_int* info );

        extern "C" void ztrtri_( const char* uplo, const char* diag, const lapack_int* n,
                lapack_complex_double* a, const lapack_int* lda, lapack_int* info );

        extern "C" void sposv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs, float* a,
                const lapack_int* lda, float* b, const lapack_int* ldb, lapack_int* info );

        extern "C" void dposv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs,
                double* a, const lapack_int* lda, double* b, const lapack_int* ldb,
                lapack_int* info );

        extern "C" void cposv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs,
                lapack_complex_float* a, const lapack_int* lda, lapack_complex_float* b,
                const lapack_int* ldb, lapack_int* info );

        extern "C" void zposv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs,
                lapack_complex_double* a, const lapack_int* lda, lapack_complex_double* b,
                const lapack_int* ldb, lapack_int* info );

        extern "C" void sgetrf_( const lapack_int* m, const lapack_int* n, float* a, const lapack_int* lda,
                lapack_int* ipiv, lapack_int* info );

        extern "C" void dgetrf_( const lapack_int* m, const lapack_int* n, double* a,
                const lapack_int* lda, lapack_int* ipiv, lapack_int* info );

        extern "C" void cgetrf_( const lapack_int* m, const lapack_int* n, lapack_complex_float* a,
                const lapack_int* lda, lapack_int* ipiv, lapack_int* info );

        extern "C" void zgetrf_( const lapack_int* m, const lapack_int* n, lapack_complex_double* a,
                const lapack_int* lda, lapack_int* ipiv, lapack_int* info );

        extern "C" void sgetri_( const lapack_int* n, float* a, const lapack_int* lda,
                const lapack_int* ipiv, float* work, const lapack_int* lwork,
                lapack_int* info );

        extern "C" void dgetri_( const lapack_int* n, double* a, const lapack_int* lda,
                const lapack_int* ipiv, double* work, const lapack_int* lwork,
                lapack_int* info );

        extern "C" void cgetri_( const lapack_int* n, lapack_complex_float* a, const lapack_int* lda,
                const lapack_int* ipiv, lapack_complex_float* work, const lapack_int* lwork,
                lapack_int* info );

        extern "C" void zgetri_( const lapack_int* n, lapack_complex_double* a, const lapack_int* lda,
                const lapack_int* ipiv, lapack_complex_double* work, const lapack_int* lwork,
                lapack_int* info );

    #endif // USE_LAPACK
#endif // USE_MKL

namespace Gadgetron
{

// following matrix computation calls MKL functions
#if defined(USE_MKL) || defined(USE_LAPACK)

bool GeneralMatrixProduct_gemm_CXFL(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B)
{
    typedef std::complex<float> T;
    try
    {
        char TA, TB;

        lapack_int lda = A.get_size(0);
        lapack_int ldb = B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.get_size(0);
        lapack_int K = A.get_size(1);

        lapack_int K2 = B.get_size(0);
        lapack_int N = B.get_size(1);

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

        GT_Complex8 alpha(1), beta(0);

        TA = 'N';
        TB = 'N';

        cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<const lapack_complex_float*>(pA), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm_CXFL(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B) ...");
        return false;
    }
    return true;
}

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoNDArray<float>& C, 
                            const hoNDArray<float>& A, bool transA, 
                            const hoNDArray<float>& B, bool transB)
{
    try
    {
        typedef float T;

        char TA, TB;

        lapack_int lda = A.get_size(0);
        lapack_int ldb = B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.get_size(0);
        lapack_int K = A.get_size(1);
        if ( transA )
        { 
            M = A.get_size(1);
            K = A.get_size(0);
        }

        lapack_int K2 = B.get_size(0);
        lapack_int N = B.get_size(1);
        if ( transB )
        {
            K2 = B.get_size(1);
            N = B.get_size(0);
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

        float alpha(1), beta(0);

        if ( transA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'T';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            sgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const float*>(pA), &lda, reinterpret_cast<const float*>(pB), &ldb, &beta, reinterpret_cast<float*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            sgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const float*>(pATmp), &lda, reinterpret_cast<const float*>(pB), &ldb, &beta, reinterpret_cast<float*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB) ...");
        return false;
    }
    return true;
}

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoNDArray<double>& C, 
                            const hoNDArray<double>& A, bool transA, 
                            const hoNDArray<double>& B, bool transB)
{
    try
    {
        typedef double T;

        char TA, TB;

        lapack_int lda = A.get_size(0);
        lapack_int ldb = B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.get_size(0);
        lapack_int K = A.get_size(1);
        if ( transA )
        { 
            M = A.get_size(1);
            K = A.get_size(0);
        }

        lapack_int K2 = B.get_size(0);
        lapack_int N = B.get_size(1);
        if ( transB )
        {
            K2 = B.get_size(1);
            N = B.get_size(0);
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

        double alpha(1), beta(0);

        if ( transA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'T';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            dgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const double*>(pA), &lda, reinterpret_cast<const double*>(pB), &ldb, &beta, reinterpret_cast<double*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            dgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const double*>(pATmp), &lda, reinterpret_cast<const double*>(pB), &ldb, &beta, reinterpret_cast<double*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB) ...");
        return false;
    }
    return true;
}

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoNDArray<GT_Complex8>& C, 
                            const hoNDArray<GT_Complex8>& A, bool transA, 
                            const hoNDArray<GT_Complex8>& B, bool transB)
{
    try
    {
        typedef GT_Complex8 T;

        char TA, TB;

        lapack_int lda = A.get_size(0);
        lapack_int ldb = B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.get_size(0);
        lapack_int K = A.get_size(1);
        if ( transA )
        { 
            M = A.get_size(1);
            K = A.get_size(0);
        }

        lapack_int K2 = B.get_size(0);
        lapack_int N = B.get_size(1);
        if ( transB )
        {
            K2 = B.get_size(1);
            N = B.get_size(0);
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

        GT_Complex8 alpha(1), beta(0);

        if ( transA )
        {
            TA = 'C';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'C';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<const lapack_complex_float*>(pA), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<lapack_complex_float*>(pATmp), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoNDArray<GT_Complex8>& C, const hoNDArray<GT_Complex8>& A, bool transA, const hoNDArray<GT_Complex8>& B, bool transB) ...");
        return false;
    }
    return true;
}

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoNDArray<GT_Complex16>& C, 
                            const hoNDArray<GT_Complex16>& A, bool transA, 
                            const hoNDArray<GT_Complex16>& B, bool transB)
{
    try
    {
        typedef GT_Complex16 T;

        char TA, TB;

        lapack_int lda = A.get_size(0);
        lapack_int ldb = B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.get_size(0);
        lapack_int K = A.get_size(1);
        if ( transA )
        { 
            M = A.get_size(1);
            K = A.get_size(0);
        }

        lapack_int K2 = B.get_size(0);
        lapack_int N = B.get_size(1);
        if ( transB )
        {
            K2 = B.get_size(1);
            N = B.get_size(0);
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

        GT_Complex16 alpha(1), beta(0);

        if ( transA )
        {
            TA = 'C';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'C';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            zgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_double*>(&alpha), reinterpret_cast<const lapack_complex_double*>(pA), &lda, reinterpret_cast<const lapack_complex_double*>(pB), &ldb, reinterpret_cast<lapack_complex_double*>(&beta), reinterpret_cast<lapack_complex_double*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            zgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_double*>(&alpha), reinterpret_cast<lapack_complex_double*>(pATmp), &lda, reinterpret_cast<const lapack_complex_double*>(pB), &ldb, reinterpret_cast<lapack_complex_double*>(&beta), reinterpret_cast<lapack_complex_double*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoNDArray<GT_Complex8>& C, const hoNDArray<GT_Complex8>& A, bool transA, const hoNDArray<GT_Complex8>& B, bool transB) ...");
        return false;
    }
    return true;
}

/// ---------------------------------------------------------------------

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoMatrix<float>& C, 
                            const hoMatrix<float>& A, bool transA, 
                            const hoMatrix<float>& B, bool transB)
{
    try
    {
        typedef float T;

        char TA, TB;

        lapack_int lda = A.rows();
        lapack_int ldb = B.rows();
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.rows();
        lapack_int K = A.cols();
        if ( transA )
        {
            M = A.cols();
            K = A.rows();
        }

        lapack_int K2 = B.rows();
        lapack_int N = B.cols();
        if ( transB )
        {
            K2 = B.cols();
            N = B.rows();
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.rows()!=M) || (C.cols()!=N) )
        {
            GADGET_CHECK_RETURN_FALSE(C.createMatrix(M, N));
        }

        T* pC = C.begin();
        lapack_int ldc = C.rows();

        float alpha(1), beta(0);

        if ( transA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'T';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            sgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const float*>(pA), &lda, reinterpret_cast<const float*>(pB), &ldb, &beta, reinterpret_cast<float*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            sgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const float*>(pATmp), &lda, reinterpret_cast<const float*>(pB), &ldb, &beta, reinterpret_cast<float*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoMatrix<float>& C, const hoMatrix<float>& A, bool transA, const hoMatrix<float>& B, bool transB) ...");
        return false;
    }
    return true;
}

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoMatrix<double>& C, 
                            const hoMatrix<double>& A, bool transA, 
                            const hoMatrix<double>& B, bool transB)
{
    try
    {
        typedef double T;

        char TA, TB;

        lapack_int lda = A.rows();
        lapack_int ldb = B.rows();
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.rows();
        lapack_int K = A.cols();
        if ( transA )
        {
            M = A.cols();
            K = A.rows();
        }

        lapack_int K2 = B.rows();
        lapack_int N = B.cols();
        if ( transB )
        {
            K2 = B.cols();
            N = B.rows();
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.rows()!=M) || (C.cols()!=N) )
        {
            GADGET_CHECK_RETURN_FALSE(C.createMatrix(M, N));
        }

        T* pC = C.begin();
        lapack_int ldc = C.rows();

        double alpha(1), beta(0);

        if ( transA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'T';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            dgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const double*>(pA), &lda, reinterpret_cast<const double*>(pB), &ldb, &beta, reinterpret_cast<double*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            dgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const double*>(pATmp), &lda, reinterpret_cast<const double*>(pB), &ldb, &beta, reinterpret_cast<double*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoMatrix<double>& C, const hoMatrix<double>& A, bool transA, const hoMatrix<double>& B, bool transB) ...");
        return false;
    }
    return true;
}

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoMatrix<GT_Complex8>& C, 
                            const hoMatrix<GT_Complex8>& A, bool transA, 
                            const hoMatrix<GT_Complex8>& B, bool transB)
{
    try
    {
        typedef GT_Complex8 T;

        char TA, TB;

        lapack_int lda = A.rows();
        lapack_int ldb = B.rows();
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.rows();
        lapack_int K = A.cols();
        if ( transA )
        {
            M = A.cols();
            K = A.rows();
        }

        lapack_int K2 = B.rows();
        lapack_int N = B.cols();
        if ( transB )
        {
            K2 = B.cols();
            N = B.rows();
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.rows()!=M) || (C.cols()!=N) )
        {
            GADGET_CHECK_RETURN_FALSE(C.createMatrix(M, N));
        }

        T* pC = C.begin();
        lapack_int ldc = C.rows();

        GT_Complex8 alpha(1), beta(0);

        if ( transA )
        {
            TA = 'C';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'C';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<const lapack_complex_float*>(pA), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<lapack_complex_float*>(pATmp), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoMatrix<GT_Complex8>& C, const hoMatrix<GT_Complex8>& A, bool transA, const hoMatrix<GT_Complex8>& B, bool transB) ...");
        return false;
    }
    return true;
}

template<> EXPORTCPUCOREMATH 
bool GeneralMatrixProduct_gemm(hoMatrix<GT_Complex16>& C, 
                            const hoMatrix<GT_Complex16>& A, bool transA, 
                            const hoMatrix<GT_Complex16>& B, bool transB)
{
    try
    {
        typedef GT_Complex16 T;

        char TA, TB;

        lapack_int lda = A.rows();
        lapack_int ldb = B.rows();
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.rows();
        lapack_int K = A.cols();
        if ( transA )
        {
            M = A.cols();
            K = A.rows();
        }

        lapack_int K2 = B.rows();
        lapack_int N = B.cols();
        if ( transB )
        {
            K2 = B.cols();
            N = B.rows();
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.rows()!=M) || (C.cols()!=N) )
        {
            GADGET_CHECK_RETURN_FALSE(C.createMatrix(M, N));
        }

        T* pC = C.begin();
        lapack_int ldc = C.rows();

        GT_Complex16 alpha(1), beta(0);

        if ( transA )
        {
            TA = 'C';
        }
        else
        {
            TA = 'N';
        }

        if ( transB )
        {
            TB = 'C';
        }
        else
        {
            TB = 'N';
        }

        if ( &A != &C )
        {
            zgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_double*>(&alpha), reinterpret_cast<const lapack_complex_double*>(pA), &lda, reinterpret_cast<const lapack_complex_double*>(pB), &ldb, reinterpret_cast<lapack_complex_double*>(&beta), reinterpret_cast<lapack_complex_double*>(pC), &ldc);
        }
        else
        {
            hoNDArray<T> aTmp(A);
            T* pATmp = aTmp.begin();
            zgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_double*>(&alpha), reinterpret_cast<lapack_complex_double*>(pATmp), &lda, reinterpret_cast<const lapack_complex_double*>(pB), &ldb, reinterpret_cast<lapack_complex_double*>(&beta), reinterpret_cast<lapack_complex_double*>(pC), &ldc);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoMatrix<GT_Complex16>& C, const hoMatrix<GT_Complex16>& A, bool transA, const hoMatrix<GT_Complex16>& B, bool transB) ...");
        return false;
    }
    return true;
}
/// ------------------------------------------------------------------------------------

template<typename T> 
bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        lapack_int info;
        lapack_int n = (lapack_int)(A.rows());
        T* pA = A.begin();
        lapack_int lda = (lapack_int)(A.rows());

        if ( typeid(T)==typeid(float) )
        {
            spotrf_(&uplo, &n, reinterpret_cast<float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dpotrf_(&uplo, &n, reinterpret_cast<double*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            cpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            zpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
        }
        else
        {
            GADGET_ERROR_MSG("CholeskyHermitianPositiveDefinite_potrf : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);

        if ( uplo == 'U' )
        {
            GADGET_CHECK_RETURN_FALSE(A.lowerTri(0));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(A.upperTri(0));
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo) ...");
        return false;
    }
    return true;
}

template<typename T> 
bool EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue)
{
    try
    {
        lapack_int M = (lapack_int)A.rows();
        GADGET_CHECK_RETURN_FALSE(A.cols() == M);

        if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
        {
            GADGET_CHECK_RETURN_FALSE(eigenValue.createMatrix(M, 1));
        }

        lapack_int info;
        char jobz = 'V';
        char uplo = 'L';
        T* pA = A.begin();
        typename realType<T>::Type* pEV = eigenValue.begin();

        //if ( typeid(T)==typeid(float) )
        //{
        //    info = LAPACKE_ssyev(LAPACK_COL_MAJOR, jobz, uplo, M, reinterpret_cast<float*>(pA), M, reinterpret_cast<float*>(pEV));
        //}
        //else if ( typeid(T)==typeid(double) )
        //{
        //    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplo, M, reinterpret_cast<double*>(pA), M, reinterpret_cast<double*>(pEV));
        //}
        //else if ( typeid(T)==typeid(GT_Complex8) )
        //{
        //    info = LAPACKE_cheev(LAPACK_COL_MAJOR, jobz, uplo, M, reinterpret_cast<lapack_complex_float*>(pA), M, reinterpret_cast<float*>(pEV));
        //}
        //else if ( typeid(T)==typeid(GT_Complex16) )
        //{
        //    info = LAPACKE_zheev(LAPACK_COL_MAJOR, jobz, uplo, M, reinterpret_cast<lapack_complex_double*>(pA), M, reinterpret_cast<double*>(pEV));
        //}
        //else
        //{
        //    GADGET_ERROR_MSG("EigenAnalysis_syev_heev : unsupported type " << typeid(T).name());
        //    return false;
        //}

        lapack_int lwork;
        lwork = M*M;

        if ( typeid(T)==typeid(float) )
        {
            hoNDArray<float> work(M, M);
            ssyev_(&jobz, &uplo, &M, reinterpret_cast<float*>(pA), &M, reinterpret_cast<float*>(pEV), work.begin(), &lwork, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            hoNDArray<double> work(M, M);
            dsyev_(&jobz, &uplo, &M, reinterpret_cast<double*>(pA), &M, reinterpret_cast<double*>(pEV), work.begin(), &lwork, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            hoNDArray<GT_Complex8> work(M, M);
            hoNDArray<float> rwork(3*M);
            cheev_(&jobz, &uplo, &M, reinterpret_cast<lapack_complex_float*>(pA), &M, reinterpret_cast<float*>(pEV), reinterpret_cast<lapack_complex_float*>(work.begin()), &lwork, rwork.begin(), &info);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            hoNDArray<GT_Complex16> work(M, M);
            hoNDArray<double> rwork(3*M);
            zheev_(&jobz, &uplo, &M, reinterpret_cast<lapack_complex_double*>(pA), &M, reinterpret_cast<double*>(pEV), reinterpret_cast<lapack_complex_double*>(work.begin()), &lwork, rwork.begin(), &info);
        }
        else
        {
            GADGET_ERROR_MSG("EigenAnalysis_syev_heev : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        lapack_int info;
        char uplo = 'L';
        lapack_int n = (lapack_int)A.rows();
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.rows();

        //if ( typeid(T)==typeid(float) )
        //{
        //    info = LAPACKE_spotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<float*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);

        //    info = LAPACKE_spotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<float*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);
        //}
        //else if ( typeid(T)==typeid(double) )
        //{
        //    info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<double*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);

        //    info = LAPACKE_dpotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<double*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);
        //}
        //else if ( typeid(T)==typeid(GT_Complex8) )
        //{
        //    info = LAPACKE_cpotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);

        //    info = LAPACKE_cpotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);
        //}
        //else if ( typeid(T)==typeid(GT_Complex16) )
        //{
        //    info = LAPACKE_zpotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_double*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);

        //    info = LAPACKE_zpotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_double*>(pA), lda);
        //    GADGET_CHECK_RETURN_FALSE(info==0);
        //}
        //else
        //{
        //    GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteInverse_potri : unsupported type " << typeid(T).name());
        //    return false;
        //}

        if ( typeid(T)==typeid(float) )
        {
            spotrf_(&uplo, &n, reinterpret_cast<float*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);

            spotri_(&uplo, &n, reinterpret_cast<float*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dpotrf_(&uplo, &n, reinterpret_cast<double*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);

            dpotri_(&uplo, &n, reinterpret_cast<double*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            cpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);

            cpotri_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            zpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);

            zpotri_(&uplo, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else
        {
            GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteInverse_potri : unsupported type " << typeid(T).name());
            return false;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A) ...");
        return false;
    }
    return true;
}

template<typename T> 
bool TriangularInverse_trtri(hoMatrix<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        lapack_int info;
        char diag = 'N';
        lapack_int n = (lapack_int)A.rows();
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.rows();

        /*if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_strtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<float*>(pA), lda);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dtrtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<double*>(pA), lda);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_ctrtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_ztrtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<lapack_complex_double*>(pA), lda);
        }
        else
        {
            GADGET_ERROR_MSG("TriangularInverse_trtri : unsupported type " << typeid(T).name());
            return false;
        }*/

        if ( typeid(T)==typeid(float) )
        {
            strtri_(&uplo, &diag, &n, reinterpret_cast<float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dtrtri_(&uplo, &diag, &n, reinterpret_cast<double*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            ctrtri_(&uplo, &diag, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            ztrtri_(&uplo, &diag, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
        }
        else
        {
            GADGET_ERROR_MSG("TriangularInverse_trtri : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in TriangularInverse_trtri(hoMatrix<float>& A, char uplo) ...");
        return false;
    }
    return true;
}

template<typename T> 
bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        if( b.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==b.rows());

        lapack_int info;
        char uplo = 'L';
        lapack_int n = (lapack_int)A.rows();
        lapack_int nrhs = (lapack_int)b.cols();
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.rows();
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.rows();

        /*if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_sposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<float*>(pA), lda, reinterpret_cast<float*>(pB), ldb);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<double*>(pA), lda, reinterpret_cast<double*>(pB), ldb);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_cposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_complex_float*>(pB), ldb);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_zposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<lapack_complex_double*>(pA), lda, reinterpret_cast<lapack_complex_double*>(pB), ldb);
        }
        else
        {
            GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteLinearSystem_posv : unsupported type " << typeid(T).name());
            return false;
        }*/

        if ( typeid(T)==typeid(float) )
        {
            sposv_(&uplo, &n, &nrhs, reinterpret_cast<float*>(pA), &lda, reinterpret_cast<float*>(pB), &ldb, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dposv_(&uplo, &n, &nrhs, reinterpret_cast<double*>(pA), &lda, reinterpret_cast<double*>(pB), &ldb, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            cposv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_complex_float*>(pB), &ldb, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            zposv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_double*>(pA), &lda, reinterpret_cast<lapack_complex_double*>(pB), &ldb, &info);
        }
        else
        {
            GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteLinearSystem_posv : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b) ...");
        return false;
    }
    return true;
}

/// Computes the LU factorization of a general m-by-n matrix
/// this function is called by general matrix inversion
template<typename T> 
bool LUFactorizationGeneralMatrix_getrf(hoMatrix<T>& A, hoNDArray<lapack_int>& ipiv)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;

        lapack_int info;
        lapack_int m = (lapack_int)A.rows();
        lapack_int n = (lapack_int)A.cols();

        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.rows();

        ipiv.create( GT_MIN(m, n) );
        lapack_int* pIPIV = ipiv.begin();

        //if ( typeid(T)==typeid(float) )
        //{
        //    info = LAPACKE_sgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else if ( typeid(T)==typeid(double) )
        //{
        //    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<double*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else if ( typeid(T)==typeid(GT_Complex8) )
        //{
        //    info = LAPACKE_cgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else if ( typeid(T)==typeid(GT_Complex16) )
        //{
        //    info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<lapack_complex_double*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else
        //{
        //    GADGET_ERROR_MSG("LUFactorizationGeneralMatrix_getrf : unsupported type " << typeid(T).name());
        //    return false;
        //}

        if ( typeid(T)==typeid(float) )
        {
            sgetrf_(&m, &n, reinterpret_cast<float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dgetrf_(&m, &n, reinterpret_cast<double*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            cgetrf_(&m, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            zgetrf_(&m, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else
        {
            GADGET_ERROR_MSG("LUFactorizationGeneralMatrix_getrf : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in LUFactorizationGeneralMatrix_getrf(hoMatrix<T>& A, hoMatrix<T>& ipiv) ...");
        return false;
    }
    return true;
}

/// Computes the inverse of an LU-factored general matrix
template<typename T> 
bool InverseGeneralMatrix_getri(hoMatrix<T>& A)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;

        lapack_int info;
        lapack_int m = (lapack_int)A.rows();
        lapack_int n = (lapack_int)A.cols();
        GADGET_CHECK_RETURN_FALSE(m==n);

        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.rows();

        hoNDArray<lapack_int> ipiv;
        GADGET_CHECK_RETURN_FALSE(LUFactorizationGeneralMatrix_getrf(A, ipiv));

        lapack_int* pIPIV = ipiv.begin();

        lapack_int lwork = m*m;

        /*if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_sgetri(LAPACK_COL_MAJOR, m, reinterpret_cast<float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dgetri(LAPACK_COL_MAJOR, m, reinterpret_cast<double*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_cgetri(LAPACK_COL_MAJOR, m, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_zgetri(LAPACK_COL_MAJOR, m, reinterpret_cast<lapack_complex_double*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        }
        else
        {
            GADGET_ERROR_MSG("InverseGeneralMatrix_getri : unsupported type " << typeid(T).name());
            return false;
        }*/

        if ( typeid(T)==typeid(float) )
        {
            hoNDArray<float> work(m, m);
            sgetri_(&m, reinterpret_cast<float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), work.begin(), &lwork, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            hoNDArray<double> work(m, m);
            dgetri_(&m, reinterpret_cast<double*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), work.begin(), &lwork, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            hoNDArray<GT_Complex8> work(m, m);
            cgetri_(&m, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), reinterpret_cast<lapack_complex_float*>(work.begin()), &lwork, &info);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            hoNDArray<GT_Complex16> work(m, m);
            zgetri_(&m, reinterpret_cast<lapack_complex_double*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), reinterpret_cast<lapack_complex_double*>(work.begin()), &lwork, &info);
        }
        else
        {
            GADGET_ERROR_MSG("InverseGeneralMatrix_getri : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in InverseGeneralMatrix_getri(hoMatrix<T>& A) ...");
        return false;
    }
    return true;
}

// ---------------------------------------------

template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<float>& A, char uplo);

template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<float>& A, hoMatrix<realType<float>::Type>& eigenValue);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<float>& A);

template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<float>& A, char uplo);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<float>& A, hoMatrix<float>& b);

template EXPORTCPUCOREMATH bool LUFactorizationGeneralMatrix_getrf(hoMatrix<float>& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<float>& A);

// ---------------------------------------------

template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<double>& A, char uplo);

template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<double>& A, hoMatrix<realType<double>::Type>& eigenValue);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<double>& A);

template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<double>& A, char uplo);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<double>& A, hoMatrix<double>& b);

template EXPORTCPUCOREMATH bool LUFactorizationGeneralMatrix_getrf(hoMatrix<double>& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<double>& A);


// ---------------------------------------------

template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<GT_Complex8>& A, char uplo);

template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<GT_Complex8>& A, hoMatrix<realType<GT_Complex8>::Type>& eigenValue);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<GT_Complex8>& A);

template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<GT_Complex8>& A, char uplo);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& b);

template EXPORTCPUCOREMATH bool LUFactorizationGeneralMatrix_getrf(hoMatrix<GT_Complex8>& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<GT_Complex8>& A);

// ---------------------------------------------

template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<GT_Complex16>& A, char uplo);

template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<GT_Complex16>& A, hoMatrix<realType<GT_Complex16>::Type>& eigenValue);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<GT_Complex16>& A);

template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<GT_Complex16>& A, char uplo);

template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& b);

template EXPORTCPUCOREMATH bool LUFactorizationGeneralMatrix_getrf(hoMatrix<GT_Complex16>& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<GT_Complex16>& A);

#else

    // matrix computation calls armadillo
    #ifdef USE_ARMADILLO

    #pragma message("Compile armadillo implementation of hoMatrix functions ... ")

    bool GeneralMatrixProduct_gemm_CXFL(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B)
    {
        return GeneralMatrixProduct_gemm(C, A, false, B, false);
    }

    template<typename T> 
    bool GeneralMatrixProduct_gemm(hoNDArray<T>& C, 
                                const hoNDArray<T>& A, bool transA, 
                                const hoNDArray<T>& B, bool transB)
    {
        try
        {
            size_t M = A.get_size(0);
            size_t K = A.get_size(1);
            if ( transA )
            { 
                M = A.get_size(1);
                K = A.get_size(0);
            }

            size_t K2 = B.get_size(0);
            size_t N = B.get_size(1);
            if ( transB )
            {
                K2 = B.get_size(1);
                N = B.get_size(0);
            }

            GADGET_CHECK_RETURN_FALSE(K==K2);
            if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
            {
                C.create(M, N);
            }

            const arma::Mat<typename stdType<T>::Type> armaA = as_arma_matrix(&A);
            const arma::Mat<typename stdType<T>::Type> armaB = as_arma_matrix(&B);
            arma::Mat<typename stdType<T>::Type> armaC = as_arma_matrix(&C);

            if ( !transA && !transB )
            {
                armaC = armaA * armaB;
            }
            else if ( transA && !transB )
            {
                armaC = arma::trans(armaA) * armaB;
            }
            else if ( !transA && transB )
            {
                armaC = armaA * arma::trans(armaB);
            }
            else
            {
                armaC = arma::trans(armaA) * arma::trans(armaB);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoNDArray<T>& C, const hoNDArray<T>& A, bool transA, const hoNDArray<T>& B, bool transB) ...");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool GeneralMatrixProduct_gemm(hoMatrix<T>& C, 
                                const hoMatrix<T>& A, bool transA, 
                                const hoMatrix<T>& B, bool transB)
    {
        try
        {
            size_t M = A.get_size(0);
            size_t K = A.get_size(1);
            if ( transA )
            { 
                M = A.get_size(1);
                K = A.get_size(0);
            }

            size_t K2 = B.get_size(0);
            size_t N = B.get_size(1);
            if ( transB )
            {
                K2 = B.get_size(1);
                N = B.get_size(0);
            }

            GADGET_CHECK_RETURN_FALSE(K==K2);
            if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
            {
                C.createMatrix(M, N);
            }

            const arma::Mat<typename stdType<T>::Type> armaA( (typename stdType<T>::Type*)(A.begin()), 
                                                               A.get_size(0), 
                                                               A.get_size(1), 
                                                               false, true);

            const arma::Mat<typename stdType<T>::Type> armaB( (typename stdType<T>::Type*)(B.begin()), 
                                                               B.get_size(0), 
                                                               B.get_size(1), 
                                                               false, true);

            arma::Mat<typename stdType<T>::Type> armaC( (typename stdType<T>::Type*)(C.begin()), 
                                                            C.get_size(0), C.get_size(1), false, true );

            if ( !transA && !transB )
            {
                armaC = armaA * armaB;
            }
            else if ( transA && !transB )
            {
                arma::Mat<typename stdType<T>::Type> AT = arma::trans(armaA);
                armaC = AT * armaB;
            }
            else if ( !transA && transB )
            {
                arma::Mat<typename stdType<T>::Type> BT = arma::trans(armaB);
                armaC = armaA * BT;
            }
            else
            {
                arma::Mat<typename stdType<T>::Type> AT = arma::trans(armaA);
                arma::Mat<typename stdType<T>::Type> BT = arma::trans(armaB);
                armaC = AT * BT;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoMatrix<T>& C, const hoMatrix<T>& A, bool transA, const hoMatrix<T>& B, bool transB) ...");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A)
    {
        try
        {
            if( A.get_number_of_elements()==0 ) return true;
            GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

            arma::Mat<typename stdType<T>::Type> armaA = as_arma_matrix(&A);
            arma::Mat<typename stdType<T>::Type> invA = arma::inv(armaA);

            memcpy(A.begin(), invA.memptr(), A.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A) ...");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo)
    {
        try
        {
            if( A.get_number_of_elements()==0 ) return true;
            GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

            arma::Mat<typename stdType<T>::Type> armaA = as_arma_matrix(&A);

            if ( uplo == 'U' )
            {
                arma::Mat<typename stdType<T>::Type> R = arma::chol(armaA);
                memcpy(A.begin(), R.begin(), A.get_number_of_bytes());
                GADGET_CHECK_RETURN_FALSE(A.lowerTri(0));
            }
            else
            {
                // need to get the lower triangle part of A
                arma::Mat<typename stdType<T>::Type> AL = arma::trimatl(armaA);
                arma::Mat<typename stdType<T>::Type> R = AL.st();

                AL = arma::chol(R);
                R = AL.st();

                memcpy(A.begin(), R.begin(), A.get_number_of_bytes());
                GADGET_CHECK_RETURN_FALSE(A.upperTri(0));
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo) ...");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue)
    {
        try
        {
            long long M = (long long)A.rows();
            GADGET_CHECK_RETURN_FALSE(A.cols() == M);

            if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
            {
                GADGET_CHECK_RETURN_FALSE(eigenValue.createMatrix(M, 1));
            }

            arma::Mat<typename stdType<T>::Type> armaA = as_arma_matrix(&A);
            arma::Col<typename realType<T>::Type > eigval = as_arma_col(&eigenValue);
            arma::Mat<typename stdType<T>::Type> eigvec;

            arma::eig_sym(eigval, eigvec, armaA);

            memcpy(A.begin(), eigvec.begin(), A.get_number_of_bytes());
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool TriangularInverse_trtri(hoMatrix<T>& A, char uplo)
    {
        try
        {
            if( A.get_number_of_elements()==0 ) return true;
            GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

            arma::Mat<typename stdType<T>::Type> armaA = as_arma_matrix(&A);

            arma::Mat<typename stdType<T>::Type> R;

            if ( uplo == 'U' )
            {
                R = arma::inv(arma::trimatu(armaA));
            }
            else
            {
                R = arma::inv(arma::trimatl(armaA));
            }

            memcpy(A.begin(), R.begin(), A.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in TriangularInverse_trtri(hoMatrix<float>& A, char uplo) ...");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b)
    {
        try
        {
            if( A.get_number_of_elements()==0 ) return true;
            if( b.get_number_of_elements()==0 ) return true;
            GADGET_CHECK_RETURN_FALSE(A.rows()==b.rows());

            arma::Mat<typename stdType<T>::Type> armaA = as_arma_matrix(&A);
            arma::Mat<typename stdType<T>::Type> armaB = as_arma_matrix(&b);

            arma::Mat<typename stdType<T>::Type> X = arma::solve(armaA, armaB);

            memcpy(b.begin(), X.begin(), b.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b) ...");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool InverseGeneralMatrix_getri(hoMatrix<T>& A)
    {
        try
        {
            if( A.get_number_of_elements()==0 ) return true;

            arma::Mat<typename stdType<T>::Type> armaA = as_arma_matrix(&A);
            arma::Mat<typename stdType<T>::Type> armaB = arma::inv(armaA);

            memcpy(A.begin(), armaB.begin(), A.get_number_of_bytes());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in InverseGeneralMatrix_getri(hoMatrix<T>& A) ...");
            return false;
        }
        return true;
    }

    // -----------------------------------------------------------------

    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB);
    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoMatrix<float>& C, const hoMatrix<float>& A, bool transA, const hoMatrix<float>& B, bool transB);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<float>& A);
    template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<float>& A, char uplo);
    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<float>& A, hoMatrix<typename realType<float>::Type>& eigenValue);
    template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<float>& A, char uplo);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<float>& A, hoMatrix<float>& b);
    template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<float>& A);

    // -----------------------------------------------------------------

    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB);
    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoMatrix<double>& C, const hoMatrix<double>& A, bool transA, const hoMatrix<double>& B, bool transB);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<double>& A);
    template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<double>& A, char uplo);
    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<double>& A, hoMatrix<typename realType<double>::Type>& eigenValue);
    template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<double>& A, char uplo);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<double>& A, hoMatrix<double>& b);
    template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<double>& A);

    // -----------------------------------------------------------------

    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoNDArray<GT_Complex8>& C, const hoNDArray<GT_Complex8>& A, bool transA, const hoNDArray<GT_Complex8>& B, bool transB);
    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoMatrix<GT_Complex8>& C, const hoMatrix<GT_Complex8>& A, bool transA, const hoMatrix<GT_Complex8>& B, bool transB);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<GT_Complex8>& A);
    template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<GT_Complex8>& A, char uplo);
    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<GT_Complex8>& A, hoMatrix<typename realType<GT_Complex8>::Type>& eigenValue);
    template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<GT_Complex8>& A, char uplo);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& b);
    template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<GT_Complex8>& A);

    // -----------------------------------------------------------------

    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoNDArray<GT_Complex16>& C, const hoNDArray<GT_Complex16>& A, bool transA, const hoNDArray<GT_Complex16>& B, bool transB);
    template EXPORTCPUCOREMATH bool GeneralMatrixProduct_gemm(hoMatrix<GT_Complex16>& C, const hoMatrix<GT_Complex16>& A, bool transA, const hoMatrix<GT_Complex16>& B, bool transB);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<GT_Complex16>& A);
    template EXPORTCPUCOREMATH bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<GT_Complex16>& A, char uplo);
    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev(hoMatrix<GT_Complex16>& A, hoMatrix<typename realType<GT_Complex16>::Type>& eigenValue);
    template EXPORTCPUCOREMATH bool TriangularInverse_trtri(hoMatrix<GT_Complex16>& A, char uplo);
    template EXPORTCPUCOREMATH bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& b);
    template EXPORTCPUCOREMATH bool InverseGeneralMatrix_getri(hoMatrix<GT_Complex16>& A);

    #endif // USE_ARMADILLO

#endif // defined(USE_MKL) || defined(USE_LAPACK)

#if defined(USE_MKL) || defined(USE_LAPACK) || defined(USE_ARMADILLO)

    template<typename T> 
    bool EigenAnalysis_syev_heev2(hoMatrix<T>& A, hoMatrix<T>& eigenValue)
    {
        try
        {
            long long M = (long long)A.rows();
            GADGET_CHECK_RETURN_FALSE(A.cols() == M);

            if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
            {
                GADGET_CHECK_RETURN_FALSE(eigenValue.createMatrix(M, 1));
            }

            hoMatrix<typename realType<T>::Type> D(M, 1);
            GADGET_CHECK_RETURN_FALSE(EigenAnalysis_syev_heev(A, D));
            //GADGET_CHECK_RETURN_FALSE(eigenValue.copyFrom(D));
            eigenValue.copyFrom(D);
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in EigenAnalysis_syev_heev2(hoMatrix<T>& A, hoMatrix<T>& eigenValue) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool SolveLinearSystem_Tikhonov(hoMatrix<T>& A, hoMatrix<T>& b, hoMatrix<T>& x, double lamda)
    {
        GADGET_CHECK_RETURN_FALSE(b.rows()==A.rows());

        hoMatrix<T> AHA(A.cols(), A.cols());
        GADGET_CHECK_RETURN_FALSE(GeneralMatrixProduct_gemm(AHA, A, true, A, false));

        GADGET_CHECK_RETURN_FALSE(x.createMatrix(A.cols(), b.cols()));
        GADGET_CHECK_RETURN_FALSE(GeneralMatrixProduct_gemm(x, A, true, b, false));

        // apply the Tikhonov regularization
        // Ideally, we shall apply the regularization is lamda*maxEigenValue
        // However, computing the maximal eigenvalue is computational intensive
        // A natural alternative is to use the trace of AHA matrix, which is the sum of all eigen values
        // Since all eigen values are positive, the lamda*maxEigenValue is only ~10-20% different from lamda*sum(all eigenValues)
        // for more information, refer to:
        // Tikhonov A.N., Goncharsky A.V., Stepanov V.V., Yagola A.G., 1995, 
        // Numerical Methods for the Solution of Ill-Posed Problems, Kluwer Academic Publishers.

        size_t col = AHA.cols();
        size_t c;

        double trA = std::abs(AHA(0, 0));
        for ( c=1; c<col; c++ )
        {
            trA += std::abs(AHA(c, c));
        }

        double value = trA*lamda/col;
        for ( c=0; c<col; c++ )
        {
            AHA(c,c) = T( (typename realType<T>::Type)(std::abs(AHA(c, c)) + value) );
        }

        GADGET_CHECK_RETURN_FALSE(SymmetricHermitianPositiveDefiniteLinearSystem_posv(AHA, x));

        return true;
    }

    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev2(hoMatrix<float>& A, hoMatrix<float>& eigenValue);
    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev2(hoMatrix<double>& A, hoMatrix<double>& eigenValue);
    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev2(hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& eigenValue);
    template EXPORTCPUCOREMATH bool EigenAnalysis_syev_heev2(hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& eigenValue);

    template EXPORTCPUCOREMATH bool SolveLinearSystem_Tikhonov(hoMatrix<float>& A, hoMatrix<float>& b, hoMatrix<float>& x, double lamda);
    template EXPORTCPUCOREMATH bool SolveLinearSystem_Tikhonov(hoMatrix<double>& A, hoMatrix<double>& b, hoMatrix<double>& x, double lamda);
    template EXPORTCPUCOREMATH bool SolveLinearSystem_Tikhonov(hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& b, hoMatrix<GT_Complex8>& x, double lamda);
    template EXPORTCPUCOREMATH bool SolveLinearSystem_Tikhonov(hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& b, hoMatrix<GT_Complex16>& x, double lamda);

#endif // defined(USE_MKL) || defined(USE_LAPACK) || defined(USE_ARMADILLO)

}
