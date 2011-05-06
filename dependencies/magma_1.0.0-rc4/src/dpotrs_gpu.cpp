/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_d
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define cublasDtrsm magmablas_dtrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_dpotrs_gpu(char uplo, magma_int_t n, magma_int_t nrhs, 
                 double *dA, magma_int_t ldda, 
                 double *dB, magma_int_t lddb, magma_int_t *info)
{
/*  -- magma (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
 
    Purpose
    =======

    DPOTRS solves a system of linear equations A*X = B with a symmetric
    positive definite matrix A using the Cholesky factorization
    A = U\*\*H*U or A = L*L\*\*H computed by DPOTRF.

    Arguments
    =========
 
    UPLO    (input) CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.

    N       (input) INTEGER
            The order of the matrix A.  N >= 0.

    NRHS    (input) INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    A       (input) DOUBLE_PRECISION array, dimension (LDA,N)
            The triangular factor U or L from the Cholesky factorization
            A = U\*\*H*U or A = L*L\*\*H, as computed by DPOTRF.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    B       (input/output) DOUBLE_PRECISION array, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    LDB     (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */

    double c_one = MAGMA_D_ONE;
    
    *info = 0 ; 
    if( (uplo != 'U') && (uplo != 'u') && (uplo != 'L') && (uplo != 'l') )
        *info = -1; 
    if( n < 0 )
        *info = -2; 
    if( nrhs < 0) 
        *info = -3; 
    if ( ldda < max(1, n) )
        *info = -5; 
    if ( lddb < max(1, n) )
        *info = -7;
    if( *info != 0 ){ 
        magma_xerbla("magma_dpotrs_gpu", info); 
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    if( (n==0) || (nrhs ==0) )
        return MAGMA_SUCCESS;

    if( (uplo=='U') || (uplo=='u') ){
        cublasDtrsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        cublasDtrsm(MagmaLeft, MagmaUpper, MagmaNoTrans,   MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
    }
    else{
        cublasDtrsm(MagmaLeft, MagmaLower, MagmaNoTrans,   MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        cublasDtrsm(MagmaLeft, MagmaLower, MagmaTrans, MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
    }

    return MAGMA_SUCCESS;
}
