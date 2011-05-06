/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated c

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_c
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define cublasCtrsm magmablas_ctrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_cgesv_gpu( magma_int_t n, magma_int_t nrhs, 
                 cuFloatComplex *dA, magma_int_t ldda,
		 magma_int_t *ipiv, 
                 cuFloatComplex *dB, magma_int_t lddb, 
                 magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    Solves a system of linear equations
      A * X = B  or  A' * X = B
    with a general N-by-N matrix A using the LU factorization computed by CGETRF_GPU.

    Arguments
    =========

    TRANS   (input) CHARACTER*1
            Specifies the form of the system of equations:
            = 'N':  A * X = B  (No transpose)
            = 'T':  A'* X = B  (Transpose)
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)

    N       (input) INTEGER
            The order of the matrix A.  N >= 0.

    NRHS    (input) INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    A       (input) COMPLEX array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by CGETRF_GPU.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    IPIV    (input) INTEGER array, dimension (N)
            The pivot indices from CGETRF; for 1<=i<=N, row i of the
            matrix was interchanged with row IPIV(i).

    B       (input/output) COMPLEX array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    LDB     (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value

    HWORK   (workspace) COMPLEX array, dimension N*NRHS
    =====================================================================    */


    cuFloatComplex c_one = MAGMA_C_ONE;
    cuFloatComplex *work = NULL;
    magma_int_t i1, i2, inc, ret;

    *info = 0;
    if (n < 0) {
        *info = -1;
    } else if (nrhs < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    } else if (lddb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return MAGMA_SUCCESS;
    }

    ret = magma_cgetrf_gpu( n, n, dA, ldda, ipiv, info);
    if ( (ret != MAGMA_SUCCESS) || ( *info != 0 ) ) {
	return ret;
    }

    work = (cuFloatComplex*)malloc(n * nrhs * sizeof(cuFloatComplex));
    if ( !work ) {
        return MAGMA_ERR_ALLOCATION;
    }
      
    /* Solve A * X = B. */
    i1  = 1;
    i2  = n;
    inc = 1;

    cublasGetMatrix( n, nrhs, sizeof(cuFloatComplex), dB, lddb, work, n);
    lapackf77_claswp(&nrhs, work, &n, &i1, &i2, ipiv, &inc);
    cublasSetMatrix( n, nrhs, sizeof(cuFloatComplex), work, n, dB, lddb);
    
    cublasCtrsm(MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit,    n, nrhs, c_one, dA, ldda, dB, lddb );
    cublasCtrsm(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );

    free(work);

    return MAGMA_SUCCESS;
}

