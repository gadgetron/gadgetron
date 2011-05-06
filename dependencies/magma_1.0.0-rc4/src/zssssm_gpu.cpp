/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @author Hatem Ltaief
       @author Mathieu Faverge

       @precisions normal z -> c d s

*/
#include "common_magma.h"

#define cublasZgemm magmablas_zgemm
//#define cublasZtrsm magmablas_ztrsm
//#define cublasZtrmm magmablas_ztrmm

extern "C" magma_int_t
magma_zssssm_gpu(char storev, magma_int_t m1, magma_int_t n1, 
                 magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                 cuDoubleComplex *dA1, magma_int_t ldda1, 
                 cuDoubleComplex *dA2, magma_int_t ldda2, 
                 cuDoubleComplex *dL1, magma_int_t lddl1, 
                 cuDoubleComplex *dL2, magma_int_t lddl2,
                 magma_int_t *IPIV, magma_int_t *info)
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

#define A1T(i,j) (dA1T + (i)*ldda1 + (j))
#define A2T(i,j) (dA2T + (i)*ldda2 + (j))
#define L1(i)    (dL1  + (i)*lddl1      )
#define L2(i,j)  (dL2  + (i)*lddl2i + (j)*lddl2j)

    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    int ip, ii, sb;
    cuDoubleComplex *dA1T, *dA2T;
    char transL;
    int lddl2i, lddl2j;

    /* Check input arguments */
    *info = 0;
    if (m1 < 0) {
        *info = -1;
    }
    else if (n1 < 0) {
        *info = -2;
    }
    else if (m2 < 0) {
        *info = -3;
    }
    else if (n2 < 0) {
        *info = -4;
    }
    else if (k < 0) {
        *info = -5;
    }
    else if (ib < 0) {
        *info = -6;
    }
    else if (ldda1 < max(1,m1)) {
        *info = -8;
    }
    else if (ldda2 < max(1,m2)) {
        *info = -10;
    }
    else if (lddl1 < max(1,ib)) {
        *info = -12;
    }
    else if (lddl2 < max(1,m2)) {
        *info = -14;
    }

    if (*info != 0)
        return MAGMA_ERR_ILLEGAL_VALUE;

    /* Quick return */
    if ((m1 == 0) || (n1 == 0) || (m2 == 0) || (n2 == 0) || (k == 0) || (ib == 0))
        return MAGMA_SUCCESS;

    if ( (storev == 'C') || (storev == 'c') ) {
        magmablas_zgetmo_in( dA1, dA1T, ldda1, m1, n1 );
        magmablas_zgetmo_in( dA2, dA2T, ldda2, m2, n2 );
        transL = MagmaTrans;
        lddl2i = 1; lddl2j = lddl2;
    } else {
        dA1T = dA1;
        dA2T = dA2;
        transL = MagmaNoTrans;
        lddl2i = lddl2; lddl2j = 1;
    }
 
    ip = 0;
    for( ii=0; ii<k; ii+=ib )
    {
        sb = min( k-ii, ib);
        
#ifndef NOSWAPBLK
        magmablas_zswapblk( 'R', n1,
                            A1T(0, 0), ldda1,
                            A2T(0, 0), ldda2,
                            ii+1, ii+ib, IPIV, 1, m1 );
#else
	{ 
	    int im;
	    for(i=0; i<ib; i++) {
		im = IPIV[ip]-1;
		
		if (im != (ii+i)) {
		    im = im - m1;
		    
		    assert( (im>=0) && (im<m1) && (im<m2) );
		    magmablas_zswap( n1, A1T(ii+i, 0), 1, A2T(im, 0), 1 );
		}
		ip++;
	    }
	}
#endif

#ifndef WITHOUTTRTRI
        /* Lower, Trans, because L1 is not transposed */
        cublasZtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                     n1, sb, 
                     c_one, L1( ii),    lddl1,
                            A1T(ii, 0), ldda1);
#else
        /* Lower, Trans, because L1 is not transposed */
        cublasZtrsm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                     n1, sb, 
                     c_one, L1( ii),    lddl1,
                            A1T(ii, 0), ldda1);
#endif

        /* Second parameter is trans because L2 is not transposed */
        cublasZgemm( MagmaNoTrans, transL, 
                     n2, m2, sb, 
                     c_neg_one, A1T(ii, 0), ldda1,
                                L2( 0, ii), lddl2, 
                     c_one,     A2T(0, 0 ), ldda2 );
    }

    if ( (storev == 'C') || (storev == 'c') ) {
        magmablas_zgetmo_out( dA1, dA1T, ldda1, m1, n1 );
        magmablas_zgetmo_out( dA2, dA2T, ldda2, m2, n2 );
    }
    return MAGMA_SUCCESS;
}

