#include "log.h"
#include "hoNDArray_linalg.h"

#include "cpp_blas.h"
#include "cpp_lapack.h"
//#include "hoNDArray_elemwise.h"
//#include "hoNDArray_reductions.h"

#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"


namespace Gadgetron
{

    template<class T>
    void gemm(hoNDArray<T>& C, const hoNDArray<T>& A, const hoNDArray<T>& B){
        gemm(C, A, false, B, false);
    }

template void gemm(hoNDArray<float>& C, const hoNDArray<float>& A, const hoNDArray<float>& B);
template void gemm(hoNDArray<double>& C, const hoNDArray<double>& A, const hoNDArray<double>& B);
template void gemm(hoNDArray<std::complex<float>>& C, const hoNDArray<std::complex<float>>& A, const hoNDArray<std::complex<float>>& B);
template void gemm(hoNDArray<std::complex<double>>& C, const hoNDArray<std::complex<double>>& A, const hoNDArray<std::complex<double>>& B);
template void gemm(hoNDArray<complext<float>>& C, const hoNDArray<complext<float>>& A, const hoNDArray<complext<float>>& B);
template void gemm(hoNDArray<complext<double>>& C, const hoNDArray<complext<double>>& A, const hoNDArray<complext<double>>& B);

template<class T> 
void gemm(hoNDArray<T>& C, const hoNDArray<T>& A, bool transA, const hoNDArray<T>& B, bool transB)

{

    GADGET_CHECK_THROW( (&C!=&A) && (&C!=&B) );

    auto lda = A.get_size(0);
    auto ldb = B.get_size(0);

    auto M = A.get_size(0);
    auto K = A.get_size(1);
    if (transA) std::swap(M,K);

    auto K2 =B.get_size(0);
    auto N = B.get_size(1);
    if (transB) std::swap(K2,N);

    GADGET_CHECK_THROW(K==K2);
    if ( (C.get_size(0)!=M) || (C.get_size(1)!=N) )
    {
        C.create(M, N);
    }
    auto ldc = C.get_size(0);
    BLAS::gemm(transA,transB,M,N,K,T(1),A.data(),lda,B.data(),ldb,T(0),C.data(),ldc);

}


template void gemm(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB);
template void gemm(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB);
template void gemm(hoNDArray<std::complex<float>>& C, const hoNDArray<std::complex<float>>& A, bool transA, const hoNDArray<std::complex<float>>& B, bool transB);
template void gemm(hoNDArray<std::complex<double>>& C, const hoNDArray<std::complex<double>>& A, bool transA, const hoNDArray<std::complex<double>>& B, bool transB);
template void gemm(hoNDArray<complext<float>>& C, const hoNDArray<complext<float>>& A, bool transA, const hoNDArray<complext<float>>& B, bool transB);
template void gemm(hoNDArray<complext<double>>& C, const hoNDArray<complext<double>>& A, bool transA, const hoNDArray<complext<double>>& B, bool transB);
/// ------------------------------------------------------------------------------------


template<class T>
void syrk(hoNDArray<T>& C, const hoNDArray<T>& A, char uplo, bool isATA)
{

        GADGET_CHECK_THROW( (&A!=&C) );
        size_t lda = (size_t)A.get_size(0);

        size_t M = (size_t)A.get_size(0);
        size_t K = (size_t)A.get_size(1);
        if ( isATA )
        { 
            M = (size_t)A.get_size(1);
            K = (size_t)A.get_size(0);
        }

        if ( (C.get_size(0)!=M) || (C.get_size(1)!=M) )
        {
            C.create(M, M);
        }

        size_t ldc = (size_t)C.get_size(0);

        BLAS::syrk((uplo == 'U'),isATA,M,K,1,A.data(),lda,0,C.data(),ldc);
}

template void syrk(hoNDArray<float>& C, const hoNDArray<float>& A, char uplo, bool isATA);
template void syrk(hoNDArray<double>& C, const hoNDArray<double>& A, char uplo, bool isATA);
template void syrk(hoNDArray<std::complex<float>>& C, const hoNDArray<std::complex<float>>& A, char uplo, bool isATA);
template void syrk(hoNDArray<std::complex<double>>& C, const hoNDArray<std::complex<double>>& A, char uplo, bool isATA);

/// ------------------------------------------------------------------------------------

template<class T>
void herk(hoNDArray< std::complex<T> >& C, const hoNDArray< std::complex<T> >& A, char uplo, bool isAHA) {

    GADGET_CHECK_THROW((&A != &C));

    size_t lda = (size_t)A.get_size(0);

    size_t N = (size_t)A.get_size(0);
    size_t K = (size_t)A.get_size(1);
    if (isAHA) {
        N = (size_t)A.get_size(1);
        K = (size_t)A.get_size(0);
    }

    if ((C.get_size(0) != N) || (C.get_size(1) != N)) {
        C.create(N, N);
    }

    size_t ldc = (size_t)C.get_size(0);

    BLAS::herk((uplo == 'U'), isAHA, N, K, 1, A.data(), lda, 0, C.data(), ldc);
}

template<class T>
void herk(hoNDArray<T>& C, const hoNDArray<T>& A, char uplo, bool isAHA)
{
    syrk(C, A, uplo, isAHA);
}

template void herk(hoNDArray<float>& C, const hoNDArray<float>& A, char uplo, bool isAHA);
template void herk(hoNDArray<double>& C, const hoNDArray<double>& A, char uplo, bool isAHA);
template void herk(hoNDArray<std::complex<float>>& C, const hoNDArray<std::complex<float>>& A, char uplo, bool isAHA);
template void herk(hoNDArray<std::complex<double>>& C, const hoNDArray<std::complex<double>>& A, char uplo, bool isAHA);

/// ------------------------------------------------------------------------------------

template <typename T>
void potrf(hoNDArray<T>& A, char uplo) {

    if (A.get_number_of_elements() == 0)
        return;
    GADGET_CHECK_THROW(A.get_size(0) == A.get_size(1));

    size_t n   = (size_t)(A.get_size(0));
    size_t lda = (size_t)(A.get_size(0));

    auto info = Lapack::potrf((uplo == 'U'), n, A.data(), lda);

    GADGET_CHECK_THROW(info == 0);

    if (uplo == 'U') {
        size_t r, c;
        for (c = 0; c < n; c++) {
            for (r = c + 1; r < n; r++) {
                A[r + c * n] = 0;
            }
        }
    } else {
        size_t r, c;
        for (r = 0; r < n; r++) {
            for (c = r + 1; c < n; c++) {
                A[r + c * n] = 0;
            }
        }
    }
}

template void potrf(hoNDArray<float>& A, char uplo);
template void potrf(hoNDArray<double>& A, char uplo);
template void potrf(hoNDArray< std::complex<float> >& A, char uplo);
template void potrf(hoNDArray< std::complex<double> >& A, char uplo);

/// ------------------------------------------------------------------------------------

template<typename T> 
void heev(hoNDArray<std::complex<T>>& A, hoNDArray<T>& eigenValue)
{
        size_t M = (size_t)A.get_size(0);
        GADGET_CHECK_THROW(A.get_size(1) == M);
        if ( (eigenValue.get_size(0)!=M) || (eigenValue.get_size(1)!=1) )
        {
            eigenValue.create(M, 1);
        }
        bool eigenvectors = true;
        bool upper = false;
        auto info = Lapack::heev(eigenvectors,upper,M,A.data(),M,eigenValue.data());
        GADGET_CHECK_THROW(info==0);
}


template<typename T>
std::enable_if_t<Core::is_floating_point_v<T>> heev(hoNDArray<T>& A, hoNDArray<T>& eigenValue)
{
    size_t M = (size_t)A.get_size(0);
    GADGET_CHECK_THROW(A.get_size(1) == M);
    if ( (eigenValue.get_size(0)!=M) || (eigenValue.get_size(1)!=1) )
    {
        eigenValue.create(M, 1);
    }
    bool eigenvectors = true;
    bool upper = false;
    auto info = Lapack::syev(eigenvectors,upper,M,A.data(),M,eigenValue.data());
    GADGET_CHECK_THROW(info==0);
}

template void heev(hoNDArray<float>& A, hoNDArray<float>& eigenValue);
template void heev(hoNDArray<double>& A, hoNDArray<double>& eigenValue);
template void heev(hoNDArray< std::complex<float> >& A, hoNDArray<float>& eigenValue);
template void heev(hoNDArray< std::complex<double> >& A, hoNDArray<double>& eigenValue);

template<typename T> 
void heev(hoNDArray< std::complex<T> >& A, hoNDArray< std::complex<T> >& eigenValue)
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

template void heev(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& eigenValue);
template void heev(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& eigenValue);

/// ------------------------------------------------------------------------------------

template<typename T> 
void potri(hoNDArray<T>& A)
{
        if( A.get_number_of_elements()==0 ) return;
        GADGET_CHECK_THROW(A.get_size(0)==A.get_size(1));
        size_t n = (size_t)A.get_size(0);
        size_t lda = (size_t)A.get_size(0);
        Lapack::potri(false,n, A.data(),lda);
}

template void potri(hoNDArray<float>& A);
template void potri(hoNDArray<double>& A);
template void potri(hoNDArray< std::complex<float> >& A);
template void potri(hoNDArray< std::complex<double> >& A);

/// ------------------------------------------------------------------------------------

template <typename T> void trtri(hoNDArray<T>& A, char uplo) {
    if (A.get_number_of_elements() == 0)
        return;
    GADGET_CHECK_THROW(A.get_size(0) == A.get_size(1));
    size_t n   = (size_t)A.get_size(0);
    size_t lda = (size_t)A.get_size(0);
    Lapack::tritri(uplo == 'U', false, n, A.data(), lda);
}

template void trtri(hoNDArray<float>& A, char uplo);
template void trtri(hoNDArray<double>& A, char uplo);
template void trtri(hoNDArray< std::complex<float> >& A, char uplo);
template void trtri(hoNDArray< std::complex<double> >& A, char uplo);

/// ------------------------------------------------------------------------------------

template<typename T>
void posv(hoNDArray<T>& A, hoNDArray<T>& b)
{
    if( A.get_number_of_elements()==0 ) return;
    if( b.get_number_of_elements()==0 ) return;
    GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

    size_t n = (size_t)A.get_size(0);
    size_t nrhs = (size_t)b.get_size(1);
    size_t lda = (size_t)A.get_size(0);
    size_t ldb = (size_t)b.get_size(0);
    /*
    We are swithcing off OpenMP threading before this call.There seems to be a bad interaction between openmp, cuda, and BLAS.
    This is a temporary fix that we should keep an eye on.
    */
//TODO Figure out what on earth is going on here. A wild guess is that someone was using a pthread OpenBlas
#ifdef USE_OMP
    int num_threads = omp_get_num_threads();
    if (!omp_in_parallel() && num_threads>1) omp_set_num_threads(1);
#endif //USE_OMP
    auto info = Lapack::posv(false,n,nrhs,A.data(),lda,b.data(),ldb);

#ifdef USE_OMP
    if (!omp_in_parallel() && num_threads>1) omp_set_num_threads(num_threads);
#endif //USE_OMP

    GADGET_CHECK_THROW(info==0);

}

template void posv(hoNDArray<float>& A, hoNDArray<float>& b);
template void posv(hoNDArray<double>& A, hoNDArray<double>& b);
template void posv(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b);
template void posv(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b);
template void posv(hoNDArray<complext<float>>& A, hoNDArray< complext<float> >& b);
template void posv(hoNDArray<complext<double>>& A, hoNDArray< complext<double> >& b);

/// ------------------------------------------------------------------------------------
namespace {

    template<class T>
    std::enable_if_t<is_complex_type_v<T>,Lapack::Int> hesv_internal(bool upper, size_t n, size_t nrhs, T* a, size_t lda, Lapack::Int* ipiv, T*b, size_t ldb){
        return Lapack::hesv(upper,n,nrhs,a,lda,ipiv,b,ldb);
    }
    template<class T>
    std::enable_if_t<!is_complex_type_v<T>,Lapack::Int> hesv_internal(bool upper, size_t n, size_t nrhs, T* a, size_t lda, Lapack::Int* ipiv, T*b, size_t ldb){
        return Lapack::sysv(upper,n,nrhs,a,lda,ipiv,b,ldb);
    }
}

    template<class T>
void hesv(hoNDArray< T>& A, hoNDArray< T >& b)
{
    if (A.get_number_of_elements() == 0)
        return;
    if (b.get_number_of_elements() == 0)
        return;
    GADGET_CHECK_THROW(A.get_size(0) == b.get_size(0));

    size_t n    = (size_t)A.get_size(0);
    size_t nrhs = (size_t)b.get_size(1);
    T* pA       = A.data();
    size_t lda  = (size_t)A.get_size(0);
    T* pB       = b.data();
    size_t ldb  = (size_t)b.get_size(0);

    hoNDArray<Lapack::Int> ipiv_array(n);
    Gadgetron::clear(ipiv_array);
    auto ipiv = ipiv_array.data();
    auto info    = hesv_internal(false, n, nrhs, pA, lda, ipiv, pB, ldb);
    GADGET_CHECK_THROW(info == 0);
}


template void hesv(hoNDArray<float>& A, hoNDArray<float>& b);
template void hesv(hoNDArray<double>& A, hoNDArray<double>& b);
template void hesv(hoNDArray<std::complex<float>>& A, hoNDArray<std::complex<float>>& b);
template void hesv(hoNDArray<std::complex<double>>& A, hoNDArray<std::complex<double>>& b);
template void hesv(hoNDArray<complext<float>>& A, hoNDArray<complext<float>>& b);
template void hesv(hoNDArray<complext<double>>& A, hoNDArray<complext<double>>& b);
/// ------------------------------------------------------------------------------------

template<class T> 
void gesv(hoNDArray<T>& A, hoNDArray<T>& b)
{

    if( A.get_number_of_elements()==0 ) return;
    if( b.get_number_of_elements()==0 ) return;
    GADGET_CHECK_THROW(A.get_size(0)==b.get_size(0));

    size_t n = (size_t)A.get_size(0);
    size_t nrhs = (size_t)b.get_size(1);
    T* pA = A.data();
    size_t lda = (size_t)A.get_size(0);
    T* pB = b.data();
    size_t ldb = (size_t)b.get_size(0);

    hoNDArray<Lapack::Int> work(n);
    Gadgetron::clear(work);
    auto ipiv = work.data();

    auto info = Lapack::gesv(n,nrhs,pA,lda,ipiv,pB,ldb);
    GADGET_CHECK_THROW(info==0);

}

template void gesv(hoNDArray<float>& A, hoNDArray<float>& b);
template void gesv(hoNDArray<double>& A, hoNDArray<double>& b);
template void gesv(hoNDArray<std::complex<float>>& A, hoNDArray<std::complex<float>>& b);
template void gesv(hoNDArray<std::complex<double>>& A, hoNDArray<std::complex<double>>& b);


/// ------------------------------------------------------------------------------------

/// Computes the LU factorization of a general m-by-n matrix
/// this function is called by general matrix inversion
template<typename T> 
void getrf(hoNDArray<T>& A, hoNDArray<Lapack::Int>& ipiv)
{
        if( A.get_number_of_elements()==0 ) return;

        size_t m = (size_t)A.get_size(0);
        size_t n = (size_t)A.get_size(1);

        T* pA = A.data();
        size_t lda = (size_t)A.get_size(0);

        ipiv.create( std::min(m, n) );
        auto pIPIV = ipiv.data();

        auto info = Lapack::getrf(m, n, pA, lda, pIPIV);
        GADGET_CHECK_THROW(info==0);
}

template void getrf(hoNDArray<float>& A, hoNDArray<Lapack::Int>& ipiv);
template void getrf(hoNDArray<double>& A, hoNDArray<Lapack::Int>& ipiv);
template void getrf(hoNDArray< std::complex<float> >& A, hoNDArray<Lapack::Int>& ipiv);
template void getrf(hoNDArray< std::complex<double> >& A, hoNDArray<Lapack::Int>& ipiv);

/// ------------------------------------------------------------------------------------

/// Computes the inverse of an LU-factored general matrix
template<typename T> 
void invert(hoNDArray<T>& A)
{
        if( A.get_number_of_elements()==0 ) return;

        size_t m = (size_t)A.get_size(0);
        size_t n = (size_t)A.get_size(1);
        GADGET_CHECK_THROW(m==n);

        T* pA = A.data();
        size_t lda = (size_t)A.get_size(0);

        hoNDArray<Lapack::Int> ipiv;
        getrf(A, ipiv);

        auto pIPIV = ipiv.data();

        auto info = Lapack::getri(n,pA,lda,pIPIV);
        GADGET_CHECK_THROW(info==0);


}

template void invert(hoNDArray<float>& A);
template void invert(hoNDArray<double>& A);
template void invert(hoNDArray< std::complex<float> >& A);
template void invert(hoNDArray< std::complex<double> >& A);

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

    double trA = abs(AHA(0, 0));
    for ( size_t c=1; c<col; c++ )
    {
        trA += abs( AHA(c, c) );
    }

    double value = trA*lamda/col;
    for (size_t c=0; c<col; c++ )
    {
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

template void SolveLinearSystem_Tikhonov(hoNDArray<float>& A, hoNDArray<float>& b, hoNDArray<float>& x, double lamda);
template void SolveLinearSystem_Tikhonov(hoNDArray<double>& A, hoNDArray<double>& b, hoNDArray<double>& x, double lamda);
template void SolveLinearSystem_Tikhonov(hoNDArray< std::complex<float> >& A, hoNDArray< std::complex<float> >& b, hoNDArray< std::complex<float> >& x, double lamda);
template void SolveLinearSystem_Tikhonov(hoNDArray< std::complex<double> >& A, hoNDArray< std::complex<double> >& b, hoNDArray< std::complex<double> >& x, double lamda);

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

        invert(ATA);

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

template void linFit(const hoNDArray<float>& x, const hoNDArray<float>& y, float& a, float& b);
template void linFit(const hoNDArray<double>& x, const hoNDArray<double>& y, double& a, double& b);
template void linFit(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, std::complex<float>& a, std::complex<float>& b);
template void linFit(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, std::complex<double>& a, std::complex<double>& b);


template <typename T>
void polyFit2(const hoNDArray<T>& x, const hoNDArray<T>& y, T& a, T& b, T& c)
{
    try
    {
        GADGET_CHECK_THROW(x.get_number_of_elements() == y.get_number_of_elements());

        size_t N = x.get_number_of_elements();

        hoNDArray<T> A;
        A.create(N, 3);

        size_t n;
        for (n=0; n<N; n++)
        {
            A(n, 0) = x(n)*x(n);
            A(n, 1) = x(n);
            A(n, 2) = 1;
        }

        hoNDArray<T> A2(A);

        hoNDArray<T> ATA(3, 3), ATy(3, 1);
        gemm(ATA, A, true, A2, false);
        gemm(ATy, A, true, y, false);

        invert(ATA);

        hoNDArray<T> res;
        res.create(3, 1);
        gemm(res, ATA, false, ATy, false);

        a = res(0);
        b = res(1);
        c = res(2);
    }
    catch (...)
    {
        GADGET_THROW("Errors in polyFit2(const hoNDArray<T>& x, const hoNDArray<T>& y, T& a, T& b, T& c) ... ");
    }
}

template void polyFit2(const hoNDArray<float>& x, const hoNDArray<float>& y, float& a, float& b, float& c);
template void polyFit2(const hoNDArray<double>& x, const hoNDArray<double>& y, double& a, double& b, double& c);
template void polyFit2(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, std::complex<float>& a, std::complex<float>& b, std::complex<float>& c);
template void polyFit2(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, std::complex<double>& a, std::complex<double>& b, std::complex<double>& c);

}
