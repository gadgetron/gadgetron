/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z 

       WARNING: this version has really poor performance 
       and cublas is prefered to this implementation. 

*/
#include "common_magma.h"
#define PRECISION_z

/*The version for fermi can be found in zsymv_fermi.cu */
#if defined(PRECISION_z) && (GPUSHMEM < 200)

#define magmablas_zsymv_130 magmablas_zsymv

#define thread_seg 128  // used in zsymv_130_kernel1 
#define threadSize 128  // used in zsymv_130_kernel2

__global__ void 
magmablas_zsymv_130_kernel1( magma_int_t m, cuDoubleComplex alpha,
                             cuDoubleComplex *A, magma_int_t lda,
                             cuDoubleComplex *x, magma_int_t incx,
                             cuDoubleComplex beta,
                             cuDoubleComplex *y, magma_int_t incy )
{
    cuDoubleComplex res = MAGMA_Z_ZERO;
    magma_int_t tid = blockIdx.x * thread_seg  + threadIdx.x;
    magma_int_t i;

    if(tid < m)
    {
#pragma unroll
        for (i=0; i<tid; i++)
        {
            res +=  A[tid + i*lda] * x[i];
        }
        y[tid] = beta * y[tid] + alpha * res;
    }
}	

__global__ void 
magmablas_zsymv_130_kernel2( magma_int_t m, cuDoubleComplex alpha,
                             cuDoubleComplex *A, magma_int_t lda,
                             cuDoubleComplex *x, magma_int_t incx,
                             cuDoubleComplex beta,
                             cuDoubleComplex *y, magma_int_t incy )
{
    __shared__ cuDoubleComplex sdata[threadSize];
    magma_int_t tx = threadIdx.x;
    magma_int_t i;
    
    cuDoubleComplex zero = MAGMA_Z_ZERO;
    cuDoubleComplex res  = MAGMA_Z_ZERO;

    magma_int_t m1 = ((m - blockIdx.y)/threadSize) * threadSize;

    for(i=blockIdx.y; i<(m1 + blockIdx.y); i+= threadSize)
    {
        res += cuConj(A[tx+i + lda*blockIdx.y]) * x[tx+i];
    }
    
    if(m > (m1 + blockIdx.y))
    {
        if( (tx + m1 + blockIdx.y) <  m )
        {
            res += cuConj(A[tx+m1+blockIdx.y + lda*blockIdx.y]) 
                *  x[tx+m1+blockIdx.y];
        }
        else 
        {
            res += zero;
        }
    }	

    sdata[tx] = res;
    __syncthreads();
    
    if(tx < 64) 
    {
        sdata[tx] += sdata[tx + 64];
    }
    __syncthreads();
	
    if(tx < 32) 
    {
        sdata[tx] += sdata[tx + 32];
        sdata[tx] += sdata[tx + 16];
        sdata[tx] += sdata[tx +  8];
        sdata[tx] += sdata[tx +  4];
        sdata[tx] += sdata[tx +  2];
        sdata[tx] += sdata[tx +  1];
    }
    if( tx == 0 ) 
    {
        y[blockIdx.y] = alpha * sdata[0] + y[blockIdx.y];
    }
}


/*************************************************************************

    Purpose
    =======

    magmablas_zsymv_130  performs the matrix-vector operation on tesla:

       y := alpha*A*x + beta*y,

    where alpha and beta are scalars, x and y are n element vectors and
    A is an n by n symmetric matrix.

    Arguments
    ==========

    UPLO   - CHARACTER*1.
             On entry, UPLO specifies whether the upper or lower
             triangular part of the array A is to be referenced as
             follows:

                UPLO = 'U' or 'u'   Only the upper triangular part of A
                                    is to be referenced.

                UPLO = 'L' or 'l'   Only the lower triangular part of A
                                    is to be referenced.

             Unchanged on exit.

    N      - INTEGER.
             On entry, N specifies the order of the matrix A.
             N must be at least zero.
             Unchanged on exit.

    ALPHA  - COMPLEX*16      .
             On entry, ALPHA specifies the scalar alpha.
             Unchanged on exit.

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
             Before entry with  UPLO = 'U' or 'u', the leading n by n
             upper triangular part of the array A must contain the upper
             triangular part of the hermitian matrix and the strictly
             lower triangular part of A is not referenced.
             Before entry with UPLO = 'L' or 'l', the leading n by n
             lower triangular part of the array A must contain the lower
             triangular part of the hermitian matrix and the strictly
             upper triangular part of A is not referenced.
             Note that the imaginary parts of the diagonal elements need
             not be set and are assumed to be zero.
             Unchanged on exit.

    LDA    - INTEGER.
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program. LDA must be at least
             max( 1, n ).
             Unchanged on exit.
             It is recommended that lda is multiple of 16. Otherwise
             performance would be deteriorated as the memory accesses
             would not be fully coalescent.

    X      - COMPLEX*16       array of dimension at least
             ( 1 + ( n - 1 )*abs( INCX ) ).
             Before entry, the incremented array X must contain the n
             element vector x.
             Unchanged on exit.

    INCX   - INTEGER.
             On entry, INCX specifies the increment for the elements of
             X. INCX must not be zero.
             Unchanged on exit.

    BETA   - COMPLEX*16      .
             On entry, BETA specifies the scalar beta. When BETA is
             supplied as zero then Y need not be set on input.
             Unchanged on exit.

    Y      - COMPLEX*16       array of dimension at least
             ( 1 + ( n - 1 )*abs( INCY ) ).
             Before entry, the incremented array Y must contain the n
             element vector y. On exit, Y is overwritten by the updated
             vector y.

    INCY   - INTEGER.
             On entry, INCY specifies the increment for the elements of
             Y. INCY must not be zero.
             Unchanged on exit.

*/

extern "C"
magma_int_t
magmablas_zsymv_130( char uplo, magma_int_t n,
                     cuDoubleComplex alpha, 
                     cuDoubleComplex *A, magma_int_t lda,
                     cuDoubleComplex *X, magma_int_t incx,
                     cuDoubleComplex beta,  
                     cuDoubleComplex *Y, magma_int_t incy)
{
    char      uplo_[2] = {uplo, 0};
    long int  upper    = lapackf77_lsame(uplo_, "U");

    /*
     * Test the input parameters.
     */
    if ((! upper) && (! lapackf77_lsame(uplo_, "L"))) {
        return -1;
    } else if ( n < 0 ) {
        return -2;
    } else if ( lda < max(1,n) ) {
        return -5;
    } else if ( incx == 0 ) {
        return -7;
    } else if ( incy == 0 ) {
        return -10;
    }

    /*
     * Quick return if possible.
     */
    if ( (n == 0) || ( MAGMA_Z_EQUAL(alpha, MAGMA_Z_ZERO) && MAGMA_Z_EQUAL(beta, MAGMA_Z_ONE) ) )
        return MAGMA_SUCCESS;

    magma_int_t blocks = (n-1)/thread_seg + 1;

    dim3 grid1( blocks, 1, 1);
    dim3 threads1(thread_seg, 1, 1);
    dim3 grid2( 1, n, 1);
    dim3 threads2(threadSize, 1, 1);
        
    /* TODO: Upper case is not implemented in MAGMA */
    if ( upper ) {
#if defined(PRECISION_z) || (defined PRECISION_c)
        fprintf(stderr, "%s: %s\n", __func__, "Upper case not implemented");
#else
        cublasZsymv(uplo, n, alpha, A, lda, X, incx, beta, Y, incy);
#endif
    }
    else
    {
        magmablas_zsymv_130_kernel1 <<< grid1, threads1 >>> 
            (n, alpha, A, lda, X, incx, beta, Y, incy);

        magmablas_zsymv_130_kernel2 <<< grid2, threads2 >>> 
            (n, alpha, A, lda, X, incx, beta, Y, incy);
    }

    return MAGMA_SUCCESS;    
}

#endif /* defined(PRECISION_z) && (GPUSHMEM < 200)*/
