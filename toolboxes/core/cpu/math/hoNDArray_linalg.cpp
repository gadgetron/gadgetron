#include "log.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

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

extern "C" void sgesv_( const lapack_int* n, const lapack_int* nrhs, float* a,
        const lapack_int* lda, lapack_int* ipiv, float* b, const lapack_int* ldb, lapack_int* info );

extern "C" void dgesv_( const lapack_int* n, const lapack_int* nrhs, double* a,
        const lapack_int* lda, lapack_int* ipiv, double* b, const lapack_int* ldb, lapack_int* info );

extern "C" void cgesv_( const lapack_int* n, const lapack_int* nrhs, lapack_complex_float* a,
        const lapack_int* lda, lapack_int* ipiv, lapack_complex_float* b, const lapack_int* ldb, lapack_int* info );

extern "C" void zgesv_( const lapack_int* n, const lapack_int* nrhs, lapack_complex_double* a,
        const lapack_int* lda, lapack_int* ipiv, lapack_complex_double* b, const lapack_int* ldb, lapack_int* info );

extern "C" void ssysv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs, float* a,
        const lapack_int* lda, lapack_int* ipiv, float* b, const lapack_int* ldb, float* work, lapack_int* lwork, lapack_int* info );

extern "C" void dsysv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs, double* a,
        const lapack_int* lda, lapack_int* ipiv, double* b, const lapack_int* ldb, double* work, lapack_int* lwork, lapack_int* info );

extern "C" void chesv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs, lapack_complex_float* a,
        const lapack_int* lda, lapack_int* ipiv, lapack_complex_float* b, const lapack_int* ldb, lapack_complex_float* work, lapack_int* lwork, lapack_int* info );

extern "C" void zhesv_( const char* uplo, const lapack_int* n, const lapack_int* nrhs, lapack_complex_double* a,
        const lapack_int* lda, lapack_int* ipiv, lapack_complex_double* b, const lapack_int* ldb, lapack_complex_double* work, lapack_int* lwork,  lapack_int* info );

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



void gemm(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B)
{
    typedef std::complex<float> T;
        char TA, TB;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) );

        lapack_int lda = (lapack_int)A.get_size(0);
        lapack_int ldb = (lapack_int)B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);

        lapack_int K2 = (lapack_int)B.get_size(0);
        lapack_int N = (lapack_int)B.get_size(1);

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

         std::complex<float>  alpha(1), beta(0);

        TA = 'N';
        TB = 'N';

        cgemm_(&TA, &TB, &M, &N, &K, reinterpret_cast<lapack_complex_float*>(&alpha), reinterpret_cast<const lapack_complex_float*>(pA), &lda, reinterpret_cast<const lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(&beta), reinterpret_cast<lapack_complex_float*>(pC), &ldc);

}

template<> EXPORTCPUCOREMATH 
void gemm(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB)
{
    try
    {
        typedef float T;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) && (&A!=&B) );

        char TA, TB;

        lapack_int lda = (lapack_int)A.get_size(0);
        lapack_int ldb = (lapack_int)B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( transA )
        { 
            M = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        lapack_int K2 = (lapack_int)B.get_size(0);
        lapack_int N = (lapack_int)B.get_size(1);
        if ( transB )
        {
            K2 = (lapack_int)B.get_size(1);
            N = (lapack_int)B.get_size(0);
        }

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

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

        lapack_int lda = (lapack_int)A.get_size(0);
        lapack_int ldb = (lapack_int)B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( transA )
        { 
            M = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        lapack_int K2 = (lapack_int)B.get_size(0);
        lapack_int N = (lapack_int)B.get_size(1);
        if ( transB )
        {
            K2 = (lapack_int)B.get_size(1);
            N = (lapack_int)B.get_size(0);
        }

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

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
        typedef  std::complex<float>  T;

        GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B)  );

        char TA, TB;

        lapack_int lda = (lapack_int)A.get_size(0);
        lapack_int ldb = (lapack_int)B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( transA )
        { 
            M = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        lapack_int K2 = (lapack_int)B.get_size(0);
        lapack_int N = (lapack_int)B.get_size(1);
        if ( transB )
        {
            K2 = (lapack_int)B.get_size(1);
            N = (lapack_int)B.get_size(0);
        }

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

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

template<> EXPORTCPUCOREMATH 
void gemm(hoNDArray< complext<float> >& C, const hoNDArray< complext<float> >& A, bool transA, const hoNDArray< complext<float> >& B, bool transB)
{
    try
    {
        typedef hoNDArray< std::complex<float> > ArrayType;
        gemm( reinterpret_cast<ArrayType&>(C), reinterpret_cast<const ArrayType&>(A), transA, reinterpret_cast<const ArrayType&>(B), transB );
    }
    catch(...)
    {
        GADGET_THROW("Errors in gemm(hoNDArray< complext<float> >& C, const hoNDArray< complext<float> >& A, bool transA, const hoNDArray< complext<float> >& B, bool transB) ...");
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

        lapack_int lda = (lapack_int)A.get_size(0);
        lapack_int ldb = (lapack_int)B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        lapack_int M = (lapack_int)A.get_size(0);
        lapack_int K = (lapack_int)A.get_size(1);
        if ( transA )
        { 
            M = (lapack_int)A.get_size(1);
            K = (lapack_int)A.get_size(0);
        }

        lapack_int K2 = (lapack_int)B.get_size(0);
        lapack_int N = (lapack_int)B.get_size(1);
        if ( transB )
        {
            K2 = (lapack_int)B.get_size(1);
            N = (lapack_int)B.get_size(0);
        }

        GADGET_CHECK_THROW(K==K2);
        if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
        {
            C.create(M, N);
        }

        T* pC = C.begin();
        lapack_int ldc = (lapack_int)C.get_size(0);

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

template<> EXPORTCPUCOREMATH 
void gemm(hoNDArray< complext<double> >& C, const hoNDArray< complext<double> >& A, bool transA, const hoNDArray< complext<double> >& B, bool transB)
{
    try
    {
        typedef hoNDArray< std::complex<double> > ArrayType;
        gemm( reinterpret_cast<ArrayType&>(C), reinterpret_cast<const ArrayType&>(A), transA, reinterpret_cast<const ArrayType&>(B), transB );
    }
    catch(...)
    {
        GADGET_THROW("Errors in gemm(hoNDArray< complext<double> >& C, const hoNDArray< complext<double> >& A, bool transA, const hoNDArray< complext<double> >& B, bool transB) ...");
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
void syrk(hoNDArray< complext<float> >& C, const hoNDArray< complext<float> >& A, char uplo, bool isATA)
{
    try
    {
        typedef  hoNDArray< std::complex<float> > ArrayType;
        syrk( reinterpret_cast<ArrayType&>(C), reinterpret_cast<const ArrayType&>(A), uplo, isATA);
    }
    catch(...)
    {
        GADGET_THROW("Errors in syrk(hoNDArray< complext<float> >& C, const hoNDArray< complext<float> >& A, char uplo, bool isATA) ...");
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
void syrk(hoNDArray< complext<double> >& C, const hoNDArray< complext<double> >& A, char uplo, bool isATA)
{
    try
    {
        typedef  hoNDArray< std::complex<double> > ArrayType;
        syrk( reinterpret_cast<ArrayType&>(C), reinterpret_cast<const ArrayType&>(A), uplo, isATA);
    }
    catch(...)
    {
        GADGET_THROW("Errors in syrk(hoNDArray< complext<double> >& C, const hoNDArray< complext<double> >& A, char uplo, bool isATA) ...");
    }
}

/// ------------------------------------------------------------------------------------

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
void herk(hoNDArray< complext<float> >& C, const hoNDArray< complext<float> >& A, char uplo, bool isATA)
{
    try
    {
        typedef  hoNDArray< std::complex<float> > ArrayType;
        herk( reinterpret_cast<ArrayType&>(C), reinterpret_cast<const ArrayType&>(A), uplo, isATA);
    }
    catch(...)
    {
        GADGET_THROW("Errors in herk(hoNDArray< complext<float> >& C, const hoNDArray< complext<float> >& A, char uplo, bool isATA) ...");
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

template<> EXPORTCPUCOREMATH 
void herk(hoNDArray< complext<double> >& C, const hoNDArray< complext<double> >& A, char uplo, bool isATA)
{
    try
    {
        typedef  hoNDArray< std::complex<double> > ArrayType;
        herk( reinterpret_cast<ArrayType&>(C), reinterpret_cast<const ArrayType&>(A), uplo, isATA);
    }
    catch(...)
    {
        GADGET_THROW("Errors in herk(hoNDArray< complext<double> >& C, const hoNDArray< complext<double> >& A, char uplo, bool isATA) ...");
    }
}

/// ------------------------------------------------------------------------------------

template<typename T> 
void potrf(hoNDArray<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==A.get_size(1));

        lapack_int info;
        lapack_int n = (lapack_int)(A.get_size(0));
        T* pA = A.begin();
        lapack_int lda = (lapack_int)(A.get_size(0));

        if ( typeid(T)==typeid(float) )
        {
            spotrf_(&uplo, &n, reinterpret_cast<float*>(pA), &lda, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dpotrf_(&uplo, &n, reinterpret_cast<double*>(pA), &lda, &info);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            cpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
            // GADGET_CHECK_THROW(A.lowerTri(0));

            size_t r, c;
            for (c=0; c<n; c++)
            {
                for (r=c+1; r<n; r++)
                {
                    pA[r + c*n] = 0;
                }
            }
        }
        else
        {
            // GADGET_CHECK_THROW(A.upperTri(0));

            size_t r, c;
            for (r=0; r<n; r++)
            {
                for (c=r+1; c<n; c++)
                {
                    pA[r + c*n] = 0;
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in potrf(hoNDArray<T>& A, char uplo) ...");
    }
}

template EXPORTCPUCOREMATH void potrf(hoNDArray<float>& A, char uplo);
template EXPORTCPUCOREMATH void potrf(hoNDArray<double>& A, char uplo);
template EXPORTCPUCOREMATH void potrf(hoNDArray< std::complex<float> >& A, char uplo);
template EXPORTCPUCOREMATH void potrf(hoNDArray< complext<float> >& A, char uplo);
template EXPORTCPUCOREMATH void potrf(hoNDArray< std::complex<double> >& A, char uplo);
template EXPORTCPUCOREMATH void potrf(hoNDArray< complext<double> >& A, char uplo);

/// ------------------------------------------------------------------------------------

template<typename T> 
void heev(hoNDArray<T>& A, hoNDArray<typename realType<T>::Type>& eigenValue)
{
    try
    {
        lapack_int M = (lapack_int)A.get_size(0);
        GADGET_CHECK_THROW(A.get_size(1) == M);

        if ( (eigenValue.get_size(0)!=M) || (eigenValue.get_size(1)!=1) )
        {
            eigenValue.create(M, 1);
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
        //else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        //{
        //    info = LAPACKE_cheev(LAPACK_COL_MAJOR, jobz, uplo, M, reinterpret_cast<lapack_complex_float*>(pA), M, reinterpret_cast<float*>(pEV));
        //}
        //else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            hoNDArray< std::complex<float> > work(M, M);
            hoNDArray<float> rwork(3*M);
            cheev_(&jobz, &uplo, &M, reinterpret_cast<lapack_complex_float*>(pA), &M, reinterpret_cast<float*>(pEV), reinterpret_cast<lapack_complex_float*>(work.begin()), &lwork, rwork.begin(), &info);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        GADGET_THROW("Errors in heev(hoNDArray<T>& A, hoNDArray<typename realType<T>::Type>& eigenValue) ... ");
    }
}

template EXPORTCPUCOREMATH void heev(hoNDArray<float>& A, hoNDArray<float>& eigenValue);
template EXPORTCPUCOREMATH void heev(hoNDArray<double>& A, hoNDArray<double>& eigenValue);
template EXPORTCPUCOREMATH void heev(hoNDArray< std::complex<float> >& A, hoNDArray<float>& eigenValue);
template EXPORTCPUCOREMATH void heev(hoNDArray< complext<float> >& A, hoNDArray<float>& eigenValue);
template EXPORTCPUCOREMATH void heev(hoNDArray< std::complex<double> >& A, hoNDArray<double>& eigenValue);
template EXPORTCPUCOREMATH void heev(hoNDArray< complext<double> >& A, hoNDArray<double>& eigenValue);

template<typename T> 
void heev(hoNDArray< std::complex<T> >& A, hoNDArray< std::complex<T> >& eigenValue)
{
    try
    {
        long long M = (long long)A.get_size(0);
        GADGET_CHECK_THROW(A.get_size(1) == M);

        if ( (eigenValue.get_size(0)!=M) || (eigenValue.get_size(1)!=1) )
        {
            eigenValue.create(M, 1);
        }

        hoNDArray<typename realType<T>::Type> D(M, 1);
        heev(A, D);
        eigenValue.copyFrom(D);
    }
    catch (...)
    {
        GADGET_THROW("Errors in heev(hoNDArray< std::complex<T> >& A, hoNDArray< std::complex<T> >& eigenValue) ... ");
    }
}

template EXPORTCPUCOREMATH void heev(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& eigenValue);
template EXPORTCPUCOREMATH void heev(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& eigenValue);

/// ------------------------------------------------------------------------------------

template<typename T> 
void potri(hoNDArray<T>& A)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==A.get_size(1));

        lapack_int info;
        char uplo = 'L';
        lapack_int n = (lapack_int)A.get_size(0);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);

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
        //else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        //{
        //    info = LAPACKE_cpotrf(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);

        //    info = LAPACKE_cpotri(LAPACK_COL_MAJOR, uplo, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        //    GADGET_CHECK_THROW(info==0);
        //}
        //else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            cpotrf_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);

            cpotri_(&uplo, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
            GADGET_CHECK_THROW(info==0);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        GADGET_THROW("Errors in potri(hoNDArray<T>& A) ...");
    }
}

template EXPORTCPUCOREMATH void potri(hoNDArray<float>& A);
template EXPORTCPUCOREMATH void potri(hoNDArray<double>& A);
template EXPORTCPUCOREMATH void potri(hoNDArray< std::complex<float> >& A);
template EXPORTCPUCOREMATH void potri(hoNDArray< complext<float> >& A);
template EXPORTCPUCOREMATH void potri(hoNDArray< std::complex<double> >& A);
template EXPORTCPUCOREMATH void potri(hoNDArray< complext<double> >& A);

/// ------------------------------------------------------------------------------------

template<typename T> 
void trtri(hoNDArray<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==A.get_size(1));

        lapack_int info;
        char diag = 'N';
        lapack_int n = (lapack_int)A.get_size(0);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);

        /*if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_strtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<float*>(pA), lda);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dtrtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<double*>(pA), lda);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            info = LAPACKE_ctrtri(LAPACK_COL_MAJOR, uplo, diag, n, reinterpret_cast<lapack_complex_float*>(pA), lda);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            ctrtri_(&uplo, &diag, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, &info);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        GADGET_THROW("Errors in trtri(hoNDArray<float>& A, char uplo) ...");
    }
}

template EXPORTCPUCOREMATH void trtri(hoNDArray<float>& A, char uplo);
template EXPORTCPUCOREMATH void trtri(hoNDArray<double>& A, char uplo);
template EXPORTCPUCOREMATH void trtri(hoNDArray< std::complex<float> >& A, char uplo);
template EXPORTCPUCOREMATH void trtri(hoNDArray< complext<float> >& A, char uplo);
template EXPORTCPUCOREMATH void trtri(hoNDArray< std::complex<double> >& A, char uplo);
template EXPORTCPUCOREMATH void trtri(hoNDArray< complext<double> >& A, char uplo);

/// ------------------------------------------------------------------------------------

template<typename T>
void posv(hoNDArray<T>& A, hoNDArray<T>& b)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info;
        char uplo = 'L';
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        /*if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_sposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<float*>(pA), lda, reinterpret_cast<float*>(pB), ldb);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<double*>(pA), lda, reinterpret_cast<double*>(pB), ldb);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            info = LAPACKE_cposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_complex_float*>(pB), ldb);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
        {
            info = LAPACKE_zposv(LAPACK_COL_MAJOR, uplo, n, nrhs, reinterpret_cast<lapack_complex_double*>(pA), lda, reinterpret_cast<lapack_complex_double*>(pB), ldb);
        }
        else
        {
            GADGET_THROW("posv : unsupported type ... ");
        }*/

        /*
        We are swithcing off OpenMP threading before this call.There seems to be a bad interaction between openmp, cuda, and BLAS.
        This is a temporary fix that we should keep an eye on.
        */

#ifdef USE_OMP
        int num_threads = omp_get_num_threads();
        if (!omp_in_parallel() && num_threads>1) omp_set_num_threads(1);
#endif //USE_OMP

        if ( typeid(T)==typeid(float) )
        {
            sposv_(&uplo, &n, &nrhs, reinterpret_cast<float*>(pA), &lda, reinterpret_cast<float*>(pB), &ldb, &info);
        }
        else if ( typeid(T)==typeid(double) )
        {
            dposv_(&uplo, &n, &nrhs, reinterpret_cast<double*>(pA), &lda, reinterpret_cast<double*>(pB), &ldb, &info);
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            cposv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_complex_float*>(pB), &ldb, &info);
        }
        else if ( (typeid(T)==typeid( std::complex<double> )) || (typeid(T)==typeid( complext<double> )) )
        {
            zposv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_double*>(pA), &lda, reinterpret_cast<lapack_complex_double*>(pB), &ldb, &info);
        }
        else
        {
#ifdef USE_OMP
            if (!omp_in_parallel() && num_threads>1) omp_set_num_threads(num_threads);
#endif //USE_OM
            GADGET_THROW("posv : unsupported type ... ");
        }

#ifdef USE_OMP
        if (!omp_in_parallel() && num_threads>1) omp_set_num_threads(num_threads);
#endif //USE_OMP

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in posv(hoNDArray<T>& A, hoNDArray<T>& b) ...");
    }
}

template EXPORTCPUCOREMATH void posv(hoNDArray<float>& A, hoNDArray<float>& b);
template EXPORTCPUCOREMATH void posv(hoNDArray<double>& A, hoNDArray<double>& b);
template EXPORTCPUCOREMATH void posv(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b);
template EXPORTCPUCOREMATH void posv(hoNDArray< complext<float> >& A, hoNDArray< complext<float> >& b);
template EXPORTCPUCOREMATH void posv(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b);
template EXPORTCPUCOREMATH void posv(hoNDArray< complext<double> >& A, hoNDArray< complext<double> >& b);

/// ------------------------------------------------------------------------------------

template<> EXPORTCPUCOREMATH
void hesv(hoNDArray< float >& A, hoNDArray< float >& b)
{
    typedef float T;
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        char uplo = 'L';
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        hoNDArray<lapack_int> ipiv_array(n);
        Gadgetron::clear(ipiv_array);
        lapack_int* ipiv = ipiv_array.begin();

        lapack_int lwork(n*n);
        hoNDArray<T> work_array(lwork);
        Gadgetron::clear(work_array);
        T* work = work_array.begin();

        ssysv_(&uplo, &n, &nrhs, reinterpret_cast<float*>(pA), &lda, ipiv, reinterpret_cast<float*>(pB), &ldb, reinterpret_cast<float*>(work), &lwork, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in hesv(hoNDArray< float >& A, hoNDArray< float >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void hesv(hoNDArray< double >& A, hoNDArray< double >& b)
{
    typedef double T;
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        char uplo = 'L';
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        hoNDArray<lapack_int> ipiv_array(n);
        Gadgetron::clear(ipiv_array);
        lapack_int* ipiv = ipiv_array.begin();

        lapack_int lwork(n*n);
        hoNDArray<T> work_array(lwork);
        Gadgetron::clear(work_array);
        T* work = work_array.begin();

        dsysv_(&uplo, &n, &nrhs, reinterpret_cast<double*>(pA), &lda, ipiv, reinterpret_cast<double*>(pB), &ldb, reinterpret_cast<double*>(work), &lwork, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in hesv(hoNDArray< double >& A, hoNDArray< double >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void hesv(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b)
{
    typedef std::complex<float> T;
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        char uplo = 'L';
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        hoNDArray<lapack_int> ipiv_array(n);
        Gadgetron::clear(ipiv_array);
        lapack_int* ipiv = ipiv_array.begin();

        lapack_int lwork(n*n);
        hoNDArray<T> work_array(lwork);
        Gadgetron::clear(work_array);
        T* work = work_array.begin();

        chesv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_float*>(pA), &lda, ipiv, reinterpret_cast<lapack_complex_float*>(pB), &ldb, reinterpret_cast<lapack_complex_float*>(work), &lwork, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in hesv(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void hesv(hoNDArray< complext<float> >& A, hoNDArray< complext<float> >& b)
{
    typedef hoNDArray< std::complex<float> > ArrayType;
    try
    {
        hesv( reinterpret_cast<ArrayType&>(A), reinterpret_cast<ArrayType&>(b) );
    }
    catch(...)
    {
        GADGET_THROW("Errors in hesv(hoNDArray< complext<float> >& A, hoNDArray< complext<float> >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void hesv(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b)
{
    typedef std::complex<double> T;
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        char uplo = 'L';
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        hoNDArray<lapack_int> ipiv_array(n);
        Gadgetron::clear(ipiv_array);
        lapack_int* ipiv = ipiv_array.begin();

        lapack_int lwork(n*n);
        hoNDArray<T> work_array(lwork);
        Gadgetron::clear(work_array);
        T* work = work_array.begin();

        zhesv_(&uplo, &n, &nrhs, reinterpret_cast<lapack_complex_double*>(pA), &lda, ipiv, reinterpret_cast<lapack_complex_double*>(pB), &ldb, reinterpret_cast<lapack_complex_double*>(work), &lwork, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in hesv(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void hesv(hoNDArray< complext<double> >& A, hoNDArray< complext<double> >& b)
{
    typedef hoNDArray< std::complex<double> > ArrayType;
    try
    {
        hesv( reinterpret_cast<ArrayType&>(A), reinterpret_cast<ArrayType&>(b) );
    }
    catch(...)
    {
        GADGET_THROW("Errors in hesv(hoNDArray< complext<double> >& A, hoNDArray< complext<double> >& b) ...");
    }
}

/// ------------------------------------------------------------------------------------

template<> EXPORTCPUCOREMATH
void gesv(hoNDArray<float>& A, hoNDArray<float>& b)
{
    typedef float T;

    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(1);

        hoNDArray<lapack_int> work(n);
        Gadgetron::clear(work);
        lapack_int* ipiv = work.begin();

        sgesv_(&n, &nrhs, reinterpret_cast<float*>(pA), &lda, ipiv, reinterpret_cast<float*>(pB), &ldb, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gesv(hoNDArray<float>& A, hoNDArray<float>& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void gesv(hoNDArray<double>& A, hoNDArray<double>& b)
{
    typedef double T;

    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        hoNDArray<lapack_int> work(n);
        Gadgetron::clear(work);
        lapack_int* ipiv = work.begin();

        dgesv_(&n, &nrhs, reinterpret_cast<double*>(pA), &lda, ipiv, reinterpret_cast<double*>(pB), &ldb, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gesv(hoNDArray<double>& A, hoNDArray<double>& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void gesv(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b)
{
    typedef std::complex<float> T;
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        hoNDArray<lapack_int> work(n);
        Gadgetron::clear(work);
        lapack_int* ipiv = work.begin();

        cgesv_(&n, &nrhs, reinterpret_cast<lapack_complex_float*>(pA), &lda, ipiv, reinterpret_cast<lapack_complex_float*>(pB), &ldb, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gesv(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void gesv(hoNDArray< complext<float> >& A, hoNDArray< complext<float> >& b)
{
    typedef hoNDArray< std::complex<float> > ArrayType;
    try
    {
        gesv( reinterpret_cast<ArrayType&>(A), reinterpret_cast<ArrayType&>(b) );
    }
    catch(...)
    {
        GADGET_THROW("Errors in gesv(hoNDArray< complext<float> >& A, hoNDArray< complext<float> >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void gesv(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b)
{
    typedef std::complex<double> T;
    try
    {
        if( A.get_number_of_elements()==0 ) return;
        if( b.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

        lapack_int info(0);
        lapack_int n = (lapack_int)A.get_size(0);
        lapack_int nrhs = (lapack_int)b.get_size(1);
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.get_size(0);

        hoNDArray<lapack_int> work(n);
        Gadgetron::clear(work);
        lapack_int* ipiv = work.begin();

        zgesv_(&n, &nrhs, reinterpret_cast<lapack_complex_double*>(pA), &lda, ipiv, reinterpret_cast<lapack_complex_double*>(pB), &ldb, &info);

        GADGET_CHECK_THROW(info==0);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gesv(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b) ...");
    }
}

template<> EXPORTCPUCOREMATH
void gesv(hoNDArray< complext<double> >& A, hoNDArray< complext<double> >& b)
{
    typedef hoNDArray< std::complex<double> > ArrayType;
    try
    {
        gesv( reinterpret_cast<ArrayType&>(A), reinterpret_cast<ArrayType&>(b) );
    }
    catch(...)
    {
        GADGET_THROW("Errors in gesv(hoNDArray< complext<double> >& A, hoNDArray< complext<double> >& b) ...");
    }
}

/// ------------------------------------------------------------------------------------

/// Computes the LU factorization of a general m-by-n matrix
/// this function is called by general matrix inversion
template<typename T> 
void getrf(hoNDArray<T>& A, hoNDArray<lapack_int>& ipiv)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;

        lapack_int info;
        lapack_int m = (lapack_int)A.get_size(0);
        lapack_int n = (lapack_int)A.get_size(1);

        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);

        ipiv.create( std::min(m, n) );
        lapack_int* pIPIV = ipiv.begin();

        //if ( typeid(T)==typeid(float) )
        //{
        //    info = LAPACKE_sgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else if ( typeid(T)==typeid(double) )
        //{
        //    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<double*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        //{
        //    info = LAPACKE_cgetrf(LAPACK_COL_MAJOR, m, n, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        //}
        //else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            cgetrf_(&m, &n, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), &info);
        }
        else if ( (typeid(T)==typeid( std::complex<double> )) || (typeid(T)==typeid( complext<double> )) )
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
        GADGET_THROW("Errors in getrf(hoNDArray<T>& A, hoNDArray<T>& ipiv) ...");
    }
}

template EXPORTCPUCOREMATH void getrf(hoNDArray<float>& A, hoNDArray<lapack_int>& ipiv);
template EXPORTCPUCOREMATH void getrf(hoNDArray<double>& A, hoNDArray<lapack_int>& ipiv);
template EXPORTCPUCOREMATH void getrf(hoNDArray< std::complex<float> >& A, hoNDArray<lapack_int>& ipiv);
template EXPORTCPUCOREMATH void getrf(hoNDArray< complext<float> >& A, hoNDArray<lapack_int>& ipiv);
template EXPORTCPUCOREMATH void getrf(hoNDArray< std::complex<double> >& A, hoNDArray<lapack_int>& ipiv);
template EXPORTCPUCOREMATH void getrf(hoNDArray< complext<double> >& A, hoNDArray<lapack_int>& ipiv);

/// ------------------------------------------------------------------------------------

/// Computes the inverse of an LU-factored general matrix
template<typename T> 
void getri(hoNDArray<T>& A)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return;

        lapack_int info;
        lapack_int m = (lapack_int)A.get_size(0);
        lapack_int n = (lapack_int)A.get_size(1);
        GADGET_CHECK_THROW(m==n);

        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.get_size(0);

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
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            info = LAPACKE_cgetri(LAPACK_COL_MAJOR, m, reinterpret_cast<lapack_complex_float*>(pA), lda, reinterpret_cast<lapack_int*>(pIPIV));
        }
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<double> )) )
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
        else if ( (typeid(T)==typeid( std::complex<float> )) || (typeid(T)==typeid( complext<float> )) )
        {
            hoNDArray< std::complex<float> > work(m, m);
            cgetri_(&m, reinterpret_cast<lapack_complex_float*>(pA), &lda, reinterpret_cast<lapack_int*>(pIPIV), reinterpret_cast<lapack_complex_float*>(work.begin()), &lwork, &info);
        }
        else if ( (typeid(T)==typeid( std::complex<double> )) || (typeid(T)==typeid( complext<double> )) )
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
        GADGET_THROW("Errors in getri(hoNDArray<T>& A) ...");
    }
}

template EXPORTCPUCOREMATH void getri(hoNDArray<float>& A);
template EXPORTCPUCOREMATH void getri(hoNDArray<double>& A);
template EXPORTCPUCOREMATH void getri(hoNDArray< std::complex<float> >& A);
template EXPORTCPUCOREMATH void getri(hoNDArray< complext<float> >& A);
template EXPORTCPUCOREMATH void getri(hoNDArray< std::complex<double> >& A);
template EXPORTCPUCOREMATH void getri(hoNDArray< complext<double> >& A);

/// ------------------------------------------------------------------------------------

template<typename T>
void SolveLinearSystem_Tikhonov(hoNDArray<T>& A, hoNDArray<T>& b, hoNDArray<T>& x, double lamda)
{
    GADGET_CHECK_THROW(b.get_size(0)==A.get_size(0));

    hoNDArray<T> AHA(A.get_size(1), A.get_size(1));
    Gadgetron::clear(AHA);

    // hoNDArray<T> ACopy(A);
    // GADGET_CHECK_THROW(gemm(AHA, ACopy, true, A, false));

    //GDEBUG_STREAM("SolveLinearSystem_Tikhonov - A = " << Gadgetron::norm2(A));
    //GDEBUG_STREAM("SolveLinearSystem_Tikhonov - b = " << Gadgetron::norm2(b));

    char uplo = 'L';
    bool isAHA = true;
    herk(AHA, A, uplo, isAHA);
    //GDEBUG_STREAM("SolveLinearSystem_Tikhonov - AHA = " << Gadgetron::norm2(AHA));

    x.create(A.get_size(1), b.get_size(1));
    gemm(x, A, true, b, false);
    //GDEBUG_STREAM("SolveLinearSystem_Tikhonov - x = " << Gadgetron::norm2(x));

    // apply the Tikhonov regularization
    // Ideally, we shall apply the regularization is lamda*maxEigenValue
    // However, computing the maximal eigenvalue is computational intensive
    // A natural alternative is to use the trace of AHA matrix, which is the sum of all eigen values
    // Since all eigen values are positive, the lamda*maxEigenValue is only ~10-20% different from lamda*sum(all eigenValues)
    // for more information, refer to:
    // Tikhonov A.N., Goncharsky A.V., Stepanov V.V., Yagola A.G., 1995,
    // Numerical Methods for the Solution of Ill-Posed Problems, Kluwer Academic Publishers.

    size_t col = AHA.get_size(0);
    size_t c;

    double trA = abs(AHA(0, 0));
    for ( c=1; c<col; c++ )
    {
        //const T v = AHA(c, c);
        //const typename realType<T>::Type rv = v.real();
        //const typename realType<T>::Type iv = v.imag();
        // trA += std::sqrt(rv*rv + iv*iv);
        trA += abs( AHA(c, c) );
    }
    //GDEBUG_STREAM("SolveLinearSystem_Tikhonov - trA = " << trA);

    double value = trA*lamda/col;
    for ( c=0; c<col; c++ )
    {
        //const T v = AHA(c, c);
        //const typename realType<T>::Type rv = v.real();
        //const typename realType<T>::Type iv = v.imag();

        //AHA(c,c) = T( (typename realType<T>::Type)( std::sqrt(rv*rv + iv*iv) + value ) );
        AHA(c,c) = T( (typename realType<T>::Type)( abs( AHA(c, c) ) + value ) );
    }

    // if the data is properly SNR unit scaled, the minimal eigen value of AHA will be around 4.0 (real and imag have noise sigma being ~1.0)
    if ( trA/col < 4.0 )
    {
        typename realType<T>::Type scalingFactor = (typename realType<T>::Type)(col*4.0/trA);
        GDEBUG_STREAM("SolveLinearSystem_Tikhonov - trA is too small : " << trA << " for matrix order : " << col);
        GDEBUG_STREAM("SolveLinearSystem_Tikhonov - scale the AHA and x by " << scalingFactor);
        Gadgetron::scal( scalingFactor, AHA);
        Gadgetron::scal( scalingFactor, x);
    }

    try
    {
        posv(AHA, x);
        //GDEBUG_STREAM("SolveLinearSystem_Tikhonov - solution = " << Gadgetron::norm2(x));
    }
    catch(...)
    {
        GERROR_STREAM("posv failed in SolveLinearSystem_Tikhonov(... ) ... ");
        GDEBUG_STREAM("A = " << Gadgetron::nrm2(A));
        GDEBUG_STREAM("b = " << Gadgetron::nrm2(b));
        GDEBUG_STREAM("AHA = " << Gadgetron::nrm2(AHA));
        GDEBUG_STREAM("trA = " << trA);
        GDEBUG_STREAM("x = " << Gadgetron::nrm2(x));

        gemm(x, A, true, b, false);
        GDEBUG_STREAM("SolveLinearSystem_Tikhonov - x = " << Gadgetron::nrm2(x));

        try
        {
            hesv(AHA, x);
        }
        catch(...)
        {
            GERROR_STREAM("hesv failed in SolveLinearSystem_Tikhonov(... ) ... ");

            gemm(x, A, true, b, false);
            GDEBUG_STREAM("SolveLinearSystem_Tikhonov - x = " << Gadgetron::nrm2(x));

            try
            {
                gesv(AHA, x);
            }
            catch(...)
            {
                GERROR_STREAM("gesv failed in SolveLinearSystem_Tikhonov(... ) ... ");
                throw;
            }
        }
    }
}

template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoNDArray<float>& A, hoNDArray<float>& b, hoNDArray<float>& x, double lamda);
template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoNDArray<double>& A, hoNDArray<double>& b, hoNDArray<double>& x, double lamda);
template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b, hoNDArray< std::complex<float> >& x, double lamda);
template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoNDArray< complext<float> >& A, hoNDArray< complext<float> >& b, hoNDArray< complext<float> >& x, double lamda);
template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b, hoNDArray< std::complex<double> >& x, double lamda);
template EXPORTCPUCOREMATH void SolveLinearSystem_Tikhonov(hoNDArray< complext<double> >& A, hoNDArray< complext<double> >& b, hoNDArray< complext<double> >& x, double lamda);

template <typename T>
void linFit(const hoNDArray<T>& x, const hoNDArray<T>& y, T& a, T& b)
{
    try
    {
        GADGET_CHECK_THROW(x.get_number_of_elements() == y.get_number_of_elements());

        size_t N = x.get_number_of_elements();

        hoNDArray<T> A;
        A.create(N, 2);

        size_t n;
        for (n=0; n<N; n++)
        {
            A(n, 0) = x(n);
            A(n, 1) = 1;
        }

        hoNDArray<T> A2(A);

        hoNDArray<T> ATA(2, 2), ATy(2, 1);
        gemm(ATA, A, true, A2, false);
        gemm(ATy, A, true, y, false);

        getri(ATA);

        hoNDArray<T> res;
        res.create(2, 1);
        gemm(res, ATA, false, ATy, false);

        a = res(0);
        b = res(1);

        /*arma::Col<T> vx = as_arma_col(&x);
        arma::Col<T> vy = as_arma_col(&y);

        arma::Col<T> P = arma::polyfit(vx, vy, 1);

        a = P(0);
        b = P(1);*/
    }
    catch (...)
    {
        GADGET_THROW("Errors in linFit(const hoNDArray<T>& x, const hoNDArray<T>& y, T& a, T& b) ... ");
    }
}

template EXPORTCPUCOREMATH void linFit(const hoNDArray<float>& x, const hoNDArray<float>& y, float& a, float& b);
template EXPORTCPUCOREMATH void linFit(const hoNDArray<double>& x, const hoNDArray<double>& y, double& a, double& b);
template EXPORTCPUCOREMATH void linFit(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, std::complex<float>& a, std::complex<float>& b);
template EXPORTCPUCOREMATH void linFit(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, std::complex<double>& a, std::complex<double>& b);


}
