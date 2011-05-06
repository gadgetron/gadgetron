/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated ds

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_d
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasStrsm magmablas_strsm
#endif
// === End defining what BLAS to use ======================================

extern "C" magma_int_t 
magma_dsgetrs_gpu(char trans, magma_int_t n, magma_int_t nrhs, 
                  float  *dA, magma_int_t ldda, 
		  magma_int_t     *ipiv, 
                  double *dB, magma_int_t lddb, 
                  double  *dX, magma_int_t lddx,
		  float  *dSX,
                  magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   
    DSGETRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general N-by-N matrix A using the LU factorization computed   
    by MAGMA_SGETRF_GPU. B and X are in DOUBLE PRECISION, and A in in SINGLE PRECISION. 
    This routine is used in the mixed precision iterative solver 
    magma_dsgesv.

    Arguments   
    =========   
    TRANS   (input) CHARACTER*1
            Specifies the form of the system of equations:
            = 'N':  A * X = B  (No transpose)
            = 'T':  A'* X = B  (Transpose)
            = 'C':  A'* X = B  (ugate transpose = Transpose)

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    dA      (input) SINGLE PRECISION array on the GPU, dimension (LDDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by CGETRF_GPU.   

    LDDA    (input) INTEGER   
            The leading dimension of the array dA.  LDDA >= max(1,N).   

    IPIV    (input) INTEGER array on the GPU, dimension (N)   
            The pivot indices from CGETRF_GPU; Row i of the   
            matrix was moved to row IPIV(i).

    dB      (input) DOUBLE PRECISION array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.

    LDDB    (input) INTEGER
            The leading dimension of the arrays X and B.  LDDB >= max(1,N).    

    dX      (output) DOUBLE PRECISION array on the GPU, dimension (LDDX, NRHS)
            On exit, the solution matrix dX.

    LDDX    (input) INTEGER   
            The leading dimension of the array dX, LDDX >= max(1,N).   

    dSX     (workspace) SINGLE PRECISION array on the GPU used as workspace,
            dimension (N, NRHS)

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
    =====================================================================    */

    float cone = MAGMA_S_ONE;
    char            trans_[2] = {trans, 0};
    long int    notran = lapackf77_lsame(trans_, "N");
    magma_int_t inc;

    *info = 0;
    if ( (! notran) &&
         (! lapackf77_lsame(trans_, "T")) &&
         (! lapackf77_lsame(trans_, "C")) ) {
      *info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (nrhs < 0) {
	*info = -3;
    } else if (ldda < n) {
	*info = -5;
    } else if (lddb < n) {
	*info = -8;
    } else if (lddx < n) {
	*info = -10;
    } else if (lddx != lddb) { /* TODO: remove it when dslaswp will have the correct interface */
	*info = -10;
    }
    if (*info != 0) {
	return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
	return MAGMA_SUCCESS;
    }
    
    if (notran) {  
      inc = 1;

      /* Get X by row applying interchanges to B and cast to single */
      /* 
       * TODO: clean dslaswp interface to have interface closer to zlaswp 
       */
      //magmablas_dslaswp(nrhs, dB, lddb, dSX, lddbx, 1, n, ipiv);
      magmablas_dslaswp(nrhs, dB, lddb, dSX, n, ipiv, inc);

      /* Solve L*X = B, overwriting B with SX. */
      cublasStrsm( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit, 
		   n, nrhs, cone, dA, ldda, dSX, n);
    
      /* Solve U*X = B, overwriting B with X. */
      cublasStrsm( MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, 
		   n, nrhs, cone, dA, ldda, dSX, n);

      magmablas_slag2d(n, nrhs, dSX, n, dX, lddx, info );
    } else {
      inc = -1;

      /* Cast the DOUBLE PRECISION RHS to SINGLE PRECISION */
      magmablas_dlag2s(n, nrhs, dB, lddb, dSX, n, info );

      /* Solve A' * X = B. */
      cublasStrsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
		  n, nrhs, cone, dA, ldda, dSX, n);
      cublasStrsm(MagmaLeft, MagmaLower, MagmaTrans, MagmaUnit,
		  n, nrhs, cone, dA, ldda, dSX, n);
      
      magmablas_dslaswp(nrhs, dX, lddx, dSX, n, ipiv, inc);
    }

    return MAGMA_SUCCESS;
    /* End of MAGMA_DSGETRS */
    
} /* magma_dsgetrs */
