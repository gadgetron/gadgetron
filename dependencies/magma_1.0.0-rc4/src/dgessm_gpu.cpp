/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @author Hatem Ltaief
       @author Mathieu Faverge

       @generated d

*/
#include "common_magma.h"

#define cublasDgemm magmablas_dgemm
//#define cublasDtrsm magmablas_dtrsm
//#define cublasDtrmm magmablas_dtrmm

extern "C" magma_int_t
magma_dgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                  magma_int_t *ipiv, 
                  double *dL1, magma_int_t lddl1, 
                  double *dL,  magma_int_t lddl, 
                  double *dA,  magma_int_t ldda, 
                  magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    SGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Arguments
    =========

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) REAL array on the GPU, dimension (LDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -7, internal GPU memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================    */

#define AT(i,j) (dAT + (i)*ldda + (j)      )
#define L(i,j)  (dL  + (i)      + (j)*lddl )
#define dL1(j)  (dL1            + (j)*lddl1)

    double c_one     = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;

    int i, s, sb;
    double *dAT;

    /* Check arguments */
    *info = 0;
    if (m < 0)
	*info = -1;
    else if (n < 0)
	*info = -2;
    else if (ldda < max(1,m))
	*info = -4;

    if (*info != 0)
        return MAGMA_SUCCESS;

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return MAGMA_SUCCESS;

    if ( (storev == 'C') || (storev == 'c') ) {
        magmablas_dgetmo_in( dA, dAT, ldda, m, n );
    } else {
        dAT = dA;
    }

    s = k / ib;
    for(i = 0; i < k; i += ib) {
        sb = min(ib, k-i);

        magmablas_dlaswp( n, dAT, ldda, i+1, i+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
        cublasDtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                     n, sb, 
                     c_one, dL1(i),   lddl1,
                            AT(i, 0), ldda);
#else
        cublasDtrsm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                     n, sb, 
                     c_one, L( i, i), lddl,
                            AT(i, 0), ldda);
#endif

        if ( (i+sb) < m) {
            cublasDgemm( MagmaNoTrans, MagmaTrans, 
                         n, m-(i+sb), sb, 
                         c_neg_one, AT(i,    0), ldda,
                                    L( i+sb, i), lddl, 
                         c_one,     AT(i+sb, 0), ldda );
        }
    }

    if ( (storev == 'C') || (storev == 'c') ) {
        magmablas_dgetmo_in( dA, dAT, ldda, m, n );
    }

    return MAGMA_SUCCESS;
    /* End of MAGMA_DGETRF_GPU */
}
