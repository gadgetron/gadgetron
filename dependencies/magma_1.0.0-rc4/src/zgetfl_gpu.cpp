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
#include <plasma.h>
#include <core_blas.h>
#include "common_magma.h"

#define cublasZgemm magmablas_zgemm
//#define cublasZtrsm magmablas_ztrsm
//#define cublasZtrmm magmablas_ztrmm

extern "C" magma_int_t
magma_zgetfl_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                  cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda,
                  cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
                  magma_int_t *ipiv, 
                  cuDoubleComplex *dwork, magma_int_t lddwork,
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


    dwork   dimension(LDDWORK, IB) to store a single panel.

    LDDWORK > max(1, ((MB+31)/32)*32)

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -7, internal GPU memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================    */

#define AT(i,j) (dAT + (i)*ib*ldda + (j)*ib)
#define hA(i,j) (hA  + (i)*ib + (j)*ib*ldha)
#define hL(j)   (hL  + (j)*ib*ldhl         )
#define hL2(j)  (hL2 + (j)*ib*ldhl         )
#define dL(j)   (dL  + (j)*ib*lddl         )
#define dL2(j)  (dL2 + (j)*ib*lddl         )

    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    magma_int_t iinfo;
    magma_int_t maxm, mindim;
    magma_int_t i, rows, cols, s, ii, sb;
    cuDoubleComplex *dAT;
#ifndef WITHOUTTRTRI
    cuDoubleComplex *dL2 = dL + ib;
    cuDoubleComplex *hL2 = hL + ib;
#endif

    /* Check arguments */
    *info = 0;
    if (m < 0)
	*info = -1;
    else if (n < 0)
	*info = -2;
    else if (ldda < max(1,m))
	*info = -4;

    if (*info != 0)
        return MAGMA_ERR_ILLEGAL_VALUE;

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return MAGMA_SUCCESS;

    /* Function Body */
    mindim = min(m, n);
    s      = mindim / ib;

    if ( ib >= mindim ) {
        /* Use CPU code. */
	lapackf77_zgetrf(&m, &n, hA, &ldha, ipiv, info);

#ifndef WITHOUTTRTRI
        CORE_zlacpy(PlasmaUpperLower, mindim, mindim, hA, ldha, hL2, ldhl );

        CORE_ztrtri( PlasmaLower, PlasmaUnit, mindim, hL2, ldhl, info );
        if (*info != 0 ) {
          fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
        }          

        cublasSetMatrix(mindim, mindim, sizeof(cuDoubleComplex), hL2, ldhl, dL2, lddl);
#endif
            
        if ( (storev == 'R') || (storev == 'r') ) {
            cublasSetMatrix(m, n, sizeof(cuDoubleComplex), hA, ldha, dwork, lddwork);
            magmablas_ztranspose( dA, ldda, dwork, lddwork, m, n );
        } else {
            cublasSetMatrix(m, n, sizeof(cuDoubleComplex), hA, ldha, dA, ldda);
        }
    }
    else {
	/* Use hybrid blocked code. */
	maxm = ((m + 31)/32)*32;

        if ( (storev == 'C') || (storev == 'c') ) {
            magmablas_zgetmo_in( dA, dAT, ldda, m, n );
        } else {
            dAT = dA;
        }
            
	for( i=0; i<s; i++ )
        {
            ii = i * ib;
            sb = min(ib, mindim-ii);
            cols = maxm - ii;

            if ( i>0 ){
                // download i-th panel
                magmablas_ztranspose( dwork, maxm, AT(0, i), ldda, sb, m );
                cublasGetMatrix( m, sb, sizeof(cuDoubleComplex), dwork, maxm, hA(0, i), ldha );
                
                // make sure that gpu queue is empty
                //cuCtxSynchronize();
#ifndef WITHOUTTRTRI
                cublasZtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n - (ii+sb), ib, 
                             c_one, dL2(i-1),    lddl, 
                                    AT(i-1,i+1), ldda );
#else
                cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n - (ii+sb), ib, 
                             c_one, AT(i-1,i-1), ldda, 
                                    AT(i-1,i+1), ldda );
#endif
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-(ii+sb), m-ii, ib, 
                             c_neg_one, AT(i-1,i+1), ldda, 
                                        AT(i,  i-1), ldda, 
                             c_one,     AT(i,  i+1), ldda );
            }

            // do the cpu part
            rows = m - ii;
            lapackf77_zgetrf( &rows, &sb, hA(i, i), &ldha, ipiv+ii, &iinfo);
            if ( (*info == 0) && (iinfo > 0) )
                *info = iinfo + ii;

            { 
                int j;
                int fin = ii + sb;
                for(j=ii ; j <fin; j++) {
                    ipiv[j] = ii + ipiv[j];
                }
            }
            magmablas_zlaswp( n-ii, AT(0, i), ldda, ii+1, ii+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
            CORE_zlacpy(PlasmaLower, sb, sb, hA(i, i), ldha, hL2(i), ldhl );
            
            CORE_ztrtri( PlasmaLower, PlasmaUnit, sb, hL2(i), ldhl, info );
            if (*info != 0 ) {
              fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
            }
            cublasSetMatrix(sb, sb, sizeof(cuDoubleComplex), hL2(i), ldhl, dL2(i), lddl);
#endif
            // upload i-th panel
            cublasSetMatrix( rows, sb, sizeof(cuDoubleComplex), hA(i, i), ldha, dwork, cols);
            magmablas_ztranspose( AT(i,i), ldda, dwork, cols, rows, sb);

            // do the small non-parallel computations
            if ( s > (i+1) ) {
#ifndef WITHOUTTRTRI
                cublasZtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             sb, sb, 
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             sb, sb, 
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             sb, m-(ii+sb), sb, 
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda, 
                             c_one,     AT(i+1, i+1), ldda );
            }
            else {
                /* Update of the last panel */
#ifndef WITHOUTTRTRI
                cublasZtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n-mindim, sb, 
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n-mindim, sb, 
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                /* m-(ii+sb) should be always 0 */
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-mindim, m-(ii+sb), sb,
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda, 
                             c_one,     AT(i+1, i+1), ldda );
            }
        }

        if ( (storev == 'C') || (storev == 'c') ) {
            magmablas_zgetmo_out( dA, dAT, ldda, m, n );
        }
    }
    return MAGMA_SUCCESS;
}
