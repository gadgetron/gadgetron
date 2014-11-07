
#include "hoNDArray_linalg.h"
#include "hoNDArray_elemwise.h"

#ifndef lapack_complex_float
    #define lapack_complex_float  std::complex<float> 
#endif // lapack_complex_float

#ifndef lapack_complex_double
    #define lapack_complex_double  std::complex<double> 
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

extern "C" void ssyrk_( const char* uplo, const char *trans, const lapack_int *n, const lapack_int *k, const float *alpha, const float *a, const lapack_int *lda, const float *beta, float *c, const lapack_int *ldc);
extern "C" void dsyrk_( const char* uplo, const char *trans, const lapack_int *n, const lapack_int *k, const double *alpha, const double *a, const lapack_int *lda, const double *beta, double *c, const lapack_int *ldc);
extern "C" void csyrk_( const char* uplo, const char *trans, const lapack_int *n, const lapack_int *k, const lapack_complex_float *alpha, const lapack_complex_float *a, const lapack_int *lda, const lapack_complex_float *beta, lapack_complex_float *c, const lapack_int *ldc);
extern "C" void zsyrk_( const char* uplo, const char *trans, const lapack_int *n, const lapack_int *k, const lapack_complex_double *alpha, const lapack_complex_double *a, const lapack_int *lda, const lapack_complex_double *beta, lapack_complex_double *c, const lapack_int *ldc);

extern "C" void cherk_( const char* uplo, const char *trans, const lapack_int *n, const lapack_int *k, const lapack_complex_float *alpha, const lapack_complex_float *a, const lapack_int *lda, const lapack_complex_float *beta, lapack_complex_float *c, const lapack_int *ldc);
extern "C" void zherk_( const char* uplo, const char *trans, const lapack_int *n, const lapack_int *k, const lapack_complex_double *alpha, const lapack_complex_double *a, const lapack_int *lda, const lapack_complex_double *beta, lapack_complex_double *c, const lapack_int *ldc);

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

namespace Gadgetron
{

// following matrix computation calls MKL functions
#if defined(USE_MKL) || defined(USE_LAPACK)

void gemm(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B)
{
    typedef std::complex<float> T;
    try
    {
        char TA, TB;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) && (&A!=&B) );

        lapack_int lda = A.get_size(0);
        lapack_int ldb = B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = A.get_size(0);
        lapack_int K = A.get_size(1);

        lapack_int K2 = B.get_size(0);
        lapack_int N = B.get_size(1);

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

         std::complex<float>  alpha(1), beta(0);

        TA = 'N';
        TB = 'N';

        cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<const lapack_complex_float*>(pA), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gemm(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void gemm(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB)
{
    try
    {
        typedef float T;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) && (&A!=&B) );

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

        GADGET_CHECK_THROW(K==K2);
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

        sgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const float*>(pA), &lda, reinterpret_cast<const float*>(pB), &ldb, &beta, reinterpret_cast<float*>(pC), &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gemm(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void gemm(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB)
{
    try
    {
        typedef double T;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) && (&A!=&B) );

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

        GADGET_CHECK_THROW(K==K2);
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

        dgemm_(&TA, &TB, &M, &N, &K, &alpha, reinterpret_cast<const double*>(pA), &lda, reinterpret_cast<const double*>(pB), &ldb, &beta, reinterpret_cast<double*>(pC), &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gemm(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void gemm(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, bool transA, const hoNDArray< std::complex<float> >& B, bool transB)
{
    try
    {
        typedef  std::complex<float>  T;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) && (&A!=&B) );

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

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

         std::complex<float>  alpha(1), beta(0);

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

        cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<const lapack_complex_float*>(pA), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gemm(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, bool transA, const hoNDArray< std::complex<float> >& B, bool transB) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void gemm(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, bool transA, const hoNDArray< std::complex<double> >& B, bool transB)
{
    try
    {
        typedef  std::complex<double>  T;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) && (&A!=&B) );

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

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = C.get_size(0);

         std::complex<double>  alpha(1), beta(0);

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

        zgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_double*>(&alpha), reinterpret_cast<const lapack_complex_double*>(pA), &lda, reinterpret_cast<const lapack_complex_double*>(pB), &ldb, reinterpret_cast<lapack_complex_double*>(&beta), reinterpret_cast<lapack_complex_double*>(pC), &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gemm(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, bool transA, const hoNDArray< std::complex<float> >& B, bool transB) ...");
    }
}

/// ------------------------------------------------------------------------------------

template<> EXPORTCPUCOREMATH 
void syrk(hoNDArray<float>& C, const hoNDArray<float>& A, char uplo, bool isATA)
{
    try
    {
        typedef float T;

        GADGET_CHECK_THROW( (&A!=&C) );

        char TA;

        lapack_int lda = (lapack_int)A.get_size(0);
        const T* pA = A.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( isATA )
        { 
            M = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        if ( (C.get_size(0)!=M) || (C.get_size(1)!=M) )
        {
            C.create(M, M);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

        float alpha(1), beta(0);

        if ( isATA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        ssyrk_(&uplo, &TA, &M, &K, &alpha, pA, &lda, &beta, pC, &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in ssyrk(hoNDArray<float>& C, const hoNDArray<float>& A, char uplo, bool isATA) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void syrk(hoNDArray<double>& C, const hoNDArray<double>& A, char uplo, bool isATA)
{
    try
    {
        typedef double T;

        GADGET_CHECK_THROW( (&A!=&C) );

        char TA;

        lapack_int lda = (lapack_int)A.get_size(0);
        const T* pA = A.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( isATA )
        { 
            M = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        if ( (C.get_size(0)!=M) || (C.get_size(1)!=M) )
        {
            C.create(M, M);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

        double alpha(1), beta(0);

        if ( isATA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        dsyrk_(&uplo, &TA, &M, &K, &alpha, pA, &lda, &beta, pC, &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in syrk(hoNDArray<double>& C, const hoNDArray<double>& A, char uplo, bool isATA) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void syrk(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, char uplo, bool isATA)
{
    try
    {
        typedef  std::complex<float>  T;

        GADGET_CHECK_THROW( (&A!=&C) );

        char TA;

        lapack_int lda = (lapack_int)A.get_size(0);
        const T* pA = A.begin(); 

        lapack_int N = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( isATA )
        { 
            N = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        GADGET_CHECK_THROW ( (C.get_size(0)==N) && (C.get_size(1)==N) );

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

        lapack_complex_float alpha(1), beta(0);

        if ( isATA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        csyrk_(&uplo, &TA, &N, &K, &alpha, pA, &lda, &beta, pC, &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in syrk(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, char uplo, bool isATA) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void syrk(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, char uplo, bool isATA)
{
    try
    {
        typedef  std::complex<double>  T;

        GADGET_CHECK_THROW( (&A!=&C) );

        char TA;

        lapack_int lda = (lapack_int)A.get_size(0);
        const T* pA = A.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( isATA )
        { 
            M = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        if ( (C.get_size(0)!=M) || (C.get_size(1)!=M) )
        {
            C.create(M, M);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

        lapack_complex_double alpha(1), beta(0);

        if ( isATA )
        {
            TA = 'T';
        }
        else
        {
            TA = 'N';
        }

        zsyrk_(&uplo, &TA, &M, &K, &alpha, pA, &lda, &beta, pC, &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in syrk(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, char uplo, bool isATA) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void herk(hoNDArray<float>& C, const hoNDArray<float>& A, char uplo, bool isAHA)
{
    syrk(C, A, uplo, isAHA);
}

template<> EXPORTCPUCOREMATH 
void herk(hoNDArray<double>& C, const hoNDArray<double>& A, char uplo, bool isAHA)
{
    syrk(C, A, uplo, isAHA);
}

template<> EXPORTCPUCOREMATH 
void herk(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, char uplo, bool isAHA)
{
    try
    {
        typedef  std::complex<float>  T;

        GADGET_CHECK_THROW( (&A!=&C) );

        char TA;

        lapack_int lda = (lapack_int)A.get_size(0);
        const T* pA = A.begin(); 

        lapack_int N = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( isAHA )
        { 
            N = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        if ( (C.get_size(0)!=N) || (C.get_size(1)!=N) )
        {
            C.create(N, N);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

        lapack_complex_float alpha(1), beta(0);

        if ( isAHA )
        {
            TA = 'C';
        }
        else
        {
            TA = 'N';
        }

        cherk_(&uplo, &TA, &N, &K, &alpha, pA, &lda, &beta, pC, &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in herk(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, char uplo, bool isAHA) ...");
    }
}

template<> EXPORTCPUCOREMATH 
void herk(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, char uplo, bool isAHA)
{
    try
    {
        typedef  std::complex<double>  T;

        GADGET_CHECK_THROW( (&A!=&C) );

        char TA;

        lapack_int lda = (lapack_int)A.get_size(0);
        const T* pA = A.begin(); 

        lapack_int N = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( isAHA )
        { 
            N = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        if ( (C.get_size(0)!=N) || (C.get_size(1)!=N) )
        {
            C.create(N, N);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

        lapack_complex_double alpha(1), beta(0);

        if ( isAHA )
        {
            TA = 'C';
        }
        else
        {
            TA = 'N';
        }

        zherk_(&uplo, &TA, &N, &K, &alpha, pA, &lda, &beta, pC, &ldc);
    }
    catch(...)
    {
        GADGET_THROW("Errors in herk(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, char uplo, bool isAHA) ...");
    }
}

/// ------------------------------------------------------------------------------------

template<typename T> 
void potrf(hoMatrix<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.rows()==A.cols());

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
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            cpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            zpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
        }
        else
        {
            GADGET_THROW("potrf : unsupported type ... ");
        }

        GADGET_CHECK_THROW(info==0);

        if ( uplo == 'U' )
        {
            GADGET_CHECK_THROW(A.lowerTri(0));
        }
        else
        {
            GADGET_CHECK_THROW(A.upperTri(0));
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in potrf(hoMatrix<T>& A, char uplo) ...");
    }
}

template<typename T> 
void heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue)
{
    try
    {
        lapack_int M = (lapack_int)A.rows();
        GADGET_CHECK_THROW(A.cols() == M);

        if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
        {
            GADGET_CHECK_THROW(eigenValue.createMatrix(M, 1));
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
        //else if ( typeid(T)==typeid( std::complex<float> ) )
        //{
        //    info = LAPACKE_cheev(LAPACK_COL_MAJOR, jobz, uplo, M, reinterpret_cast<lapack_complex_float*>(pA), M, reinterpret_cast<float*>(pEV));
        //}
        //else if ( typeid(T)==typeid( std::complex<double> ) )
        //{
        //    info = LAPACKE_zheev(LAPACK_COL_MAJOR, jobz, uplo, M, reinterpret_cast<lapack_complex_double*>(pA), M, reinterpret_cast<double*>(pEV));
        //}
        //else
        //{
        //    GADGET_THROW("heev : unsupported type " << typeid(T).name());
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
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            hoNDArray< std::complex<float> > work(M, M);
            hoNDArray<float> rwork(3*M);
            cheev_(&jobz, &uplo, &M, reinterpret_cast<lapack_complex_float*>(pA), &M, reinterpret_cast<float*>(pEV), reinterpret_cast<lapack_complex_float*>(work.begin()), &lwork, rwork.begin(), &info);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            hoNDArray< std::complex<double> > work(M, M);
            hoNDArray<double> rwork(3*M);
            zheev_(&jobz, &uplo, &M, reinterpret_cast<lapack_complex_double*>(pA), &M, reinterpret_cast<double*>(pEV), reinterpret_cast<lapack_complex_double*>(work.begin()), &lwork, rwork.begin(), &info);
        }
        else
        {
            GADGET_THROW("heev : unsupported type ... ");
        }

        GADGET_CHECK_THROW(info==0);
    }
    catch (...)
    {
        GADGET_THROW("Errors in heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue) ... ");
    }
}

template<typename T> 
void potri(hoMatrix<T>& A)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.rows()==A.cols());

        lapack_int info;
        char uplo = 'L';
        lapack_int n = (lapack_int)A.rows();
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.rows();

        //if ( typeid(T)==typeid(float) )
        //{
        //    info = LAPACKE_spotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<float*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);

        //    info = LAPACKE_spotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<float*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);
        //}
        //else if ( typeid(T)==typeid(double) )
        //{
        //    info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<double*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);

        //    info = LAPACKE_dpotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<double*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);
        //}
        //else if ( typeid(T)==typeid( std::complex<float> ) )
        //{
        //    info = LAPACKE_cpotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);

        //    info = LAPACKE_cpotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);
        //}
        //else if ( typeid(T)==typeid( std::complex<double> ) )
        //{
        //    info = LAPACKE_zpotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_double*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);

        //    info = LAPACKE_zpotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_double*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);
        //}
        //else
        //{
        //    GADGET_THROW("potri : unsupported type " << typeid(T).name());
        //}

        if ( typeid(T)==typeid(float) )
        {
            spotrf_(&uplo, &n, reinterpret_cast<float*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);

            spotri_(&uplo, &n, reinterpret_cast<float*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dpotrf_(&uplo, &n, reinterpret_cast<double*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);

            dpotri_(&uplo, &n, reinterpret_cast<double*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);
        }
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            cpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);

            cpotri_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            zpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);

            zpotri_(&uplo, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);
        }
        else
        {
            GADGET_THROW("potri : unsupported type ... ");
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in potri(hoMatrix<T>& A) ...");
    }
}

template<typename T> 
void trtri(hoMatrix<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.rows()==A.cols());

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
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            info = LAPACKE_ctrtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            info = LAPACKE_ztrtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<lapack_complex_double*>(pA), lda);
        }
        else
        {
            GADGET_THROW("trtri : unsupported type " << typeid(T).name());
        }*/

        if ( typeid(T)==typeid(float) )
        {
            strtri_(&uplo, &diag, &n, reinterpret_cast<float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dtrtri_(&uplo, &diag, &n, reinterpret_cast<double*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            ctrtri_(&uplo, &diag, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            ztrtri_(&uplo, &diag, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, &info);
        }
        else
        {
            GADGET_THROW("trtri : unsupported type ... ");
        }

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in trtri(hoMatrix<float>& A, char uplo) ...");
    }
}

template<typename T> 
void posv(hoMatrix<T>& A, hoMatrix<T>& b)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.rows()==b.rows());

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
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            info = LAPACKE_cposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_complex_float*>(pB), ldb);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            info = LAPACKE_zposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<lapack_complex_double*>(pA), lda, reinterpret_cast<lapack_complex_double*>(pB), ldb);
        }
        else
        {
            GADGET_THROW("posv : unsupported type " << typeid(T).name());
        }*/

        if ( typeid(T)==typeid(float) )
        {
            sposv_(&uplo, &n, &nrhs, reinterpret_cast<float*>(pA), &lda, reinterpret_cast<float*>(pB), &ldb, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dposv_(&uplo, &n, &nrhs, reinterpret_cast<double*>(pA), &lda, reinterpret_cast<double*>(pB), &ldb, &info);
        }
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            cposv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_complex_float*>(pB), &ldb, &info);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            zposv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_double*>(pA), &lda, reinterpret_cast<lapack_complex_double*>(pB), &ldb, &info);
        }
        else
        {
            GADGET_THROW("posv : unsupported type ... ");
        }

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in posv(hoMatrix<T>& A, hoMatrix<T>& b) ...");
    }
}

/// Computes the LU factorization of a general m-by-n matrix
/// this function is called by general matrix inversion
template<typename T> 
void getrf(hoMatrix<T>& A, hoNDArray<lapack_int>& ipiv)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;

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
        //else if ( typeid(T)==typeid( std::complex<float> ) )
        //{
        //    info = LAPACKE_cgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else if ( typeid(T)==typeid( std::complex<double> ) )
        //{
        //    info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<lapack_complex_double*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else
        //{
        //    GADGET_THROW("getrf : unsupported type " << typeid(T).name());
        //}

        if ( typeid(T)==typeid(float) )
        {
            sgetrf_(&m, &n, reinterpret_cast<float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dgetrf_(&m, &n, reinterpret_cast<double*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            cgetrf_(&m, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            zgetrf_(&m, &n, reinterpret_cast<lapack_complex_double*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else
        {
            GADGET_THROW("getrf : unsupported type ... ");
        }

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in getrf(hoMatrix<T>& A, hoMatrix<T>& ipiv) ...");
    }
}

/// Computes the inverse of an LU-factored general matrix
template<typename T> 
void getri(hoMatrix<T>& A)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;

        lapack_int info;
        lapack_int m = (lapack_int)A.rows();
        lapack_int n = (lapack_int)A.cols();
        GADGET_CHECK_THROW(m==n);

        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.rows();

        hoNDArray<lapack_int> ipiv;
        getrf(A, ipiv);

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
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            info = LAPACKE_cgetri(LAPACK_COL_MAJOR, m, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            info = LAPACKE_zgetri(LAPACK_COL_MAJOR, m, reinterpret_cast<lapack_complex_double*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        }
        else
        {
            GADGET_THROW("getri : unsupported type " << typeid(T).name());
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
        else if ( typeid(T)==typeid( std::complex<float> ) )
        {
            hoNDArray< std::complex<float> > work(m, m);
            cgetri_(&m, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), reinterpret_cast<lapack_complex_float*>(work.begin()), &lwork, &info);
        }
        else if ( typeid(T)==typeid( std::complex<double> ) )
        {
            hoNDArray< std::complex<double> > work(m, m);
            zgetri_(&m, reinterpret_cast<lapack_complex_double*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), reinterpret_cast<lapack_complex_double*>(work.begin()), &lwork, &info);
        }
        else
        {
            GADGET_THROW("getri : unsupported type ... ");
        }

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in getri(hoMatrix<T>& A) ...");
    }
}

// ---------------------------------------------

template EXPORTCPUCOREMATH void potrf(hoMatrix<float>& A, char uplo);

template EXPORTCPUCOREMATH void heev(hoMatrix<float>& A, hoMatrix<float>& eigenValue);

template EXPORTCPUCOREMATH void potri(hoMatrix<float>& A);

template EXPORTCPUCOREMATH void trtri(hoMatrix<float>& A, char uplo);

template EXPORTCPUCOREMATH void posv(hoMatrix<float>& A, hoMatrix<float>& b);

template EXPORTCPUCOREMATH void getrf(hoMatrix<float>& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH void getri(hoMatrix<float>& A);

// ---------------------------------------------

template EXPORTCPUCOREMATH void potrf(hoMatrix<double>& A, char uplo);

template EXPORTCPUCOREMATH void heev(hoMatrix<double>& A, hoMatrix<double>& eigenValue);

template EXPORTCPUCOREMATH void potri(hoMatrix<double>& A);

template EXPORTCPUCOREMATH void trtri(hoMatrix<double>& A, char uplo);

template EXPORTCPUCOREMATH void posv(hoMatrix<double>& A, hoMatrix<double>& b);

template EXPORTCPUCOREMATH void getrf(hoMatrix<double>& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH void getri(hoMatrix<double>& A);


// ---------------------------------------------

template EXPORTCPUCOREMATH void potrf(hoMatrix< std::complex<float> >& A, char uplo);

template EXPORTCPUCOREMATH void heev(hoMatrix< std::complex<float> >& A, hoMatrix<float>& eigenValue);

template EXPORTCPUCOREMATH void potri(hoMatrix< std::complex<float> >& A);

template EXPORTCPUCOREMATH void trtri(hoMatrix< std::complex<float> >& A, char uplo);

template EXPORTCPUCOREMATH void posv(hoMatrix< std::complex<float> >& A, hoMatrix< std::complex<float> >& b);

template EXPORTCPUCOREMATH void getrf(hoMatrix< std::complex<float> >& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH void getri(hoMatrix< std::complex<float> >& A);

// ---------------------------------------------

template EXPORTCPUCOREMATH void potrf(hoMatrix< std::complex<double> >& A, char uplo);

template EXPORTCPUCOREMATH void heev(hoMatrix< std::complex<double> >& A, hoMatrix<double>& eigenValue);

template EXPORTCPUCOREMATH void potri(hoMatrix< std::complex<double> >& A);

template EXPORTCPUCOREMATH void trtri(hoMatrix< std::complex<double> >& A, char uplo);

template EXPORTCPUCOREMATH void posv(hoMatrix< std::complex<double> >& A, hoMatrix< std::complex<double> >& b);

template EXPORTCPUCOREMATH void getrf(hoMatrix< std::complex<double> >& A, hoNDArray<lapack_int>& ipiv);

template EXPORTCPUCOREMATH void getri(hoMatrix< std::complex<double> >& A);

#endif // defined(USE_MKL) || defined(USE_LAPACK)

#if defined(USE_MKL) || defined(USE_LAPACK)

    template<typename T> 
    void heev(hoMatrix< std::complex<T> >& A, hoMatrix< std::complex<T> >& eigenValue)
    {
        try
        {
            long long M = (long long)A.rows();
            GADGET_CHECK_THROW(A.cols() == M);

            if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
            {
                GADGET_CHECK_THROW(eigenValue.createMatrix(M, 1));
            }

            hoMatrix<typename realType<T>::Type> D(M, 1);
            heev(A, D);
            eigenValue.copyFrom(D);
        }
        catch (...)
        {
            GADGET_THROW("Errors in heev(hoMatrix< std::complex<T> >& A, hoMatrix< std::complex<T> >& eigenValue) ... ");
        }
    }

    template<typename T> 
    void SolveLinearSystem_Tikhonov(hoMatrix<T>& A, hoMatrix<T>& b, hoMatrix<T>& x, double lamda)
    {
        GADGET_CHECK_THROW(b.rows()==A.rows());

        hoMatrix<T> AHA(A.cols(), A.cols());
        Gadgetron::clear(AHA);

        // hoMatrix<T> ACopy(A);
        // GADGET_CHECK_THROW(gemm(AHA, ACopy, true, A, false));

        char uplo = 'L';
        bool isAHA = true;
        herk(AHA, A, uplo, isAHA);

        GADGET_CHECK_THROW(x.createMatrix(A.cols(), b.cols()));
        gemm(x, A, true, b, false);

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

        posv(AHA, x);
    }

    template EXPORTCPUCOREMATH void heev(hoMatrix< std::complex<float> >& A, hoMatrix< std::complex<float> >& eigenValue);
    template EXPORTCPUCOREMATH void heev(hoMatrix< std::complex<double> >& A, hoMatrix< std::complex<double> >& eigenValue);

    template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoMatrix<float>& A, hoMatrix<float>& b, hoMatrix<float>& x, double lamda);
    template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoMatrix<double>& A, hoMatrix<double>& b, hoMatrix<double>& x, double lamda);
    template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoMatrix< std::complex<float> >& A, hoMatrix< std::complex<float> >& b, hoMatrix< std::complex<float> >& x, double lamda);
    template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoMatrix< std::complex<double> >& A, hoMatrix< std::complex<double> >& b, hoMatrix< std::complex<double> >& x, double lamda);

#endif // defined(USE_MKL) || defined(USE_LAPACK)

}
