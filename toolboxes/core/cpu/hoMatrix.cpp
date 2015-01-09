
#include "hoMatrix.h"

namespace Gadgetron
{

// C = A*B
bool GeneralMatrixProduct(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB)
{
    try
    {
        typedef float T;

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

        const T* pA = A.begin();
        const T* pB = B.begin();
        T* pC = C.begin();

        size_t m, n, k;

        if ( !transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*pB[k+n*K];
                    }
                }
            }
        }

        if ( transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[k+m*K]*pB[k+n*K];
                    }
                }
            }
        }

        if ( !transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*pB[n+k*K];
                    }
                }
            }
        }

        if ( transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[k+m*K]*pB[n+k*K];
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GeneralMatrixProduct(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB) ...");
        return false;
    }
    return true;
}

bool GeneralMatrixProduct(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB)
{
    try
    {
        typedef double T;

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

        const T* pA = A.begin();
        const T* pB = B.begin();
        T* pC = C.begin();

        size_t m, n, k;

        if ( !transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*pB[k+n*K];
                    }
                }
            }
        }

        if ( transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[k+m*K]*pB[k+n*K];
                    }
                }
            }
        }

        if ( !transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*pB[n+k*K];
                    }
                }
            }
        }

        if ( transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[k+m*K]*pB[n+k*K];
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GeneralMatrixProduct(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB) ...");
        return false;
    }
    return true;
}

bool GeneralMatrixProduct(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, bool transA, const hoNDArray< std::complex<float> >& B, bool transB)
{
    try
    {
        typedef  std::complex<float>  T;

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

        const T* pA = A.begin();
        const T* pB = B.begin();
        T* pC = C.begin();

        size_t m, n, k;

        if ( !transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*pB[k+n*K];
                    }
                }
            }
        }

        if ( transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += std::conj(pA[k+m*K])*pB[k+n*K];
                    }
                }
            }
        }

        if ( !transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*std::conj(pB[n+k*K]);
                    }
                }
            }
        }

        if ( transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += std::conj(pA[k+m*K])*std::conj(pB[n+k*K]);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GeneralMatrixProduct(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, bool transA, const hoNDArray< std::complex<float> >& B, bool transB) ...");
        return false;
    }
    return true;
}

bool GeneralMatrixProduct(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, bool transA, const hoNDArray< std::complex<double> >& B, bool transB)
{
    try
    {
        typedef  std::complex<double>  T;

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

        const T* pA = A.begin();
        const T* pB = B.begin();
        T* pC = C.begin();

        size_t m, n, k;

        if ( !transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*pB[k+n*K];
                    }
                }
            }
        }

        if ( transA && !transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += std::conj(pA[k+m*K])*pB[k+n*K];
                    }
                }
            }
        }

        if ( !transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += pA[m+k*M]*std::conj(pB[n+k*K]);
                    }
                }
            }
        }

        if ( transA && transB )
        {
            for ( m=0; m<M; m++ )
            {
                for ( n=0; n<N; n++ )
                {
                    pC[m+n*M] = 0;
                    for ( k=0; k<K; k++ )
                    {
                        pC[m+n*M] += std::conj(pA[k+m*K])*std::conj(pB[n+k*K]);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GeneralMatrixProduct(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, bool transA, const hoNDArray< std::complex<double> >& B, bool transB) ...");
        return false;
    }
    return true;
}

}
