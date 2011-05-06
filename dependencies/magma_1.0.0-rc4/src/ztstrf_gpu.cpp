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
#include "common_magma.h"

#define cublasZgemm magmablas_zgemm
//#define cublasZtrsm magmablas_ztrsm
//#define cublasZtrmm magmablas_ztrmm

extern "C" magma_int_t
magma_ztstrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                  cuDoubleComplex *hU, magma_int_t ldhu, cuDoubleComplex *dU, magma_int_t lddu, 
                  cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda, 
                  cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
                  magma_int_t *ipiv, 
                  cuDoubleComplex *hwork, magma_int_t ldhwork, cuDoubleComplex *dwork, magma_int_t lddwork,
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

#define UT(i,j) (dUT + (i)*ib*lddu + (j)*ib )
#define AT(i,j) (dAT + (i)*ib*ldda + (j)*ib )
#define L(i)    (dL  + (i)*ib*lddl          )
#define L2(i)   (dL2 + (i)*ib*lddl          )
#define hU(i,j) (hU  + (j)*ib*ldhu + (i)*ib )
#define hA(i,j) (hA  + (j)*ib*ldha + (i)*ib )
#define hL(i)   (hL  + (i)*ib*ldhl          )
#define hL2(i)  (hL2 + (i)*ib*ldhl          )

    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    int iinfo = 0;
    int maxm, mindim;
    int i, j, im, s, ip, ii, p = 1;
    cuDoubleComplex *dAT, *dUT;
    cuDoubleComplex *dAp, *dUp;
#ifndef WITHOUTTRTRI
    cuDoubleComplex *dL2 = dL + ib;
    cuDoubleComplex *hL2 = hL + ib;
    p = 2;
#endif

    /* Check input arguments */
    *info = 0;
    if (m < 0) {
        *info = -1;
    }
    else if (n < 0) {
        *info = -2;
    }
    else if (ib < 0) {
        *info = -3;
    }
    else if ((lddu < max(1,m)) && (m > 0)) {
        *info = -6;
    }
    else if ((ldda < max(1,m)) && (m > 0)) {
        *info = -8;
    }
    else if ((lddl < max(1,ib)) && (ib > 0)) {
        *info = -10;
    }

    if (*info != 0)
        return MAGMA_ERR_ILLEGAL_VALUE;

    /* quick return */
    if ((m == 0) || (n == 0) || (ib == 0))
        return MAGMA_SUCCESS;

    ip = 0;

    /* Function Body */
    mindim = min(m, n);
    s      = mindim / ib;

    if ( ib >= mindim ) {
	/* Use CPU code. */
        CORE_ztstrf(m, n, ib, nb,
                    hU, ldhu,
                    hA, ldha,
                    hL, ldhl,
                    ipiv,
                    hwork, ldhwork,
                    info);

#ifndef WITHOUTTRTRI
        CORE_zlacpy( PlasmaUpperLower, mindim, mindim, hL, ldhl, hL2, ldhl );
        CORE_ztrtri( PlasmaLower, PlasmaUnit, mindim, hL2, ldhl, info );
        if (*info != 0 ) {
          fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
        }          
#endif

	if ( (storev == 'R') || (storev == 'r') ) {
            cublasSetMatrix(m, n, sizeof(cuDoubleComplex), hU, ldhu, dwork, lddwork);
            magmablas_ztranspose( dU, lddu, dwork, lddwork, m, n );

            cublasSetMatrix(m, n, sizeof(cuDoubleComplex), hA, ldha, dwork, lddwork);
            magmablas_ztranspose( dA, ldda, dwork, lddwork, m, n );
        } else {
            cublasSetMatrix(m, n, sizeof(cuDoubleComplex), hU, ldhu, dU, lddu);
            cublasSetMatrix(m, n, sizeof(cuDoubleComplex), hA, ldha, dA, ldda);
        }
	cublasSetMatrix(p*ib, n, sizeof(cuDoubleComplex), hL, ldhl, dL, lddl);
            
    }
    else {
	/* Use hybrid blocked code. */
	maxm = ((m + 31)/32)*32;

        if ( (storev == 'C') || (storev == 'c') ) {
            magmablas_zgetmo_in( dU, dUT, lddu, m,  n );
            magmablas_zgetmo_in( dA, dAT, ldda, m,  n );
        } else {
            dUT = dU; dAT = dA;
        }
        dAp = dwork;
        dUp = dAp + ib*lddwork;

        ip = 0;
	for( i=0; i<s; i++ )
        {
            ii = i * ib;
            sb = min(mindim-ii, ib);
            
            if ( i>0 ){
                // download i-th panel
                magmablas_ztranspose( dUp, lddu, UT(0, i), lddu, sb, ii );
                magmablas_ztranspose( dAp, ldda, AT(0, i), ldda, sb, m  );
                
                cublasGetMatrix( ii, sb, sizeof(cuDoubleComplex), dUp, lddu, hU(0, i), ldhu);
                cublasGetMatrix( m,  sb, sizeof(cuDoubleComplex), dAp, ldda, hA(0, i), ldha);
                
                // make sure that gpu queue is empty
                //cuCtxSynchronize();
                
#ifndef WITHOUTTRTRI
                cublasZtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n-(ii+sb), ib, 
                             c_one, L2(i-1),      lddl,
                                    UT(i-1, i+1), lddu);
#else
                cublasZtrsm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n-(ii+sb), ib, 
                             c_one, L(i-1),       lddl,
                                    UT(i-1, i+1), lddu);
#endif
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-(ii+sb), m, ib,
                             c_neg_one, UT(i-1, i+1), lddu, 
                                        AT(0,   i-1), ldda,
                             c_one,     AT(0,   i+1), ldda );
            }

            // do the cpu part
            CORE_ztstrf(m, sb, ib, nb,
                        hU(i, i), ldhu, 
                        hA(0, i), ldha,
                        hL(i),    ldhl, 
                        ipiv+ii, 
                        hwork, ldhwork, 
                        info);

            if ( (*info == 0) && (iinfo > 0) )
                *info = iinfo + ii;
            
            // Need to swap betw U and A
#ifndef NOSWAPBLK
            magmablas_zswapblk( 'R', n-(ii+sb),
                                UT(i, i+1), lddu,
                                AT(0, i+1), ldda,
                                1, sb, ipiv+ii, 1, nb );

            for(j=0; j<ib; j++) {
                im = ipiv[ip]-1;
                if ( im == j ) {
                    ipiv[ip] += ii;
                }
                ip++;
            }
#else
            for(j=0; j<ib; j++) {
                im = ipiv[ip]-1;
                if ( im != (j) ) {
                    im = im - nb;
                    assert( (im>=0) && (im<m) );
                    magmablas_zswap( n-(ii+sb), UT(i, i+1)+j*lddu, 1, AT(0, i+1)+im*ldda, 1 );
                } else {
                    ipiv[ip] += ii;
                }
                ip++;
            }
#endif

#ifndef WITHOUTTRTRI
	    CORE_zlacpy( PlasmaUpperLower, sb, sb, hL(i), ldhl, hL2(i), ldhl );
            CORE_ztrtri( PlasmaLower, PlasmaUnit, sb, hL2(i), ldhl, info );
            if (*info != 0 ) {
              fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
            }   
#endif
            // upload i-th panel
            cublasSetMatrix( sb,   sb, sizeof(cuDoubleComplex), hU(i, i), ldhu, dUp,  lddu);
            cublasSetMatrix( m,    sb, sizeof(cuDoubleComplex), hA(0, i), ldha, dAp,  ldda);
            cublasSetMatrix( p*ib, sb, sizeof(cuDoubleComplex), hL(i),    ldhl, L(i), lddl);
            magmablas_ztranspose( UT(i, i), lddu, dUp, lddu, sb, sb);
            magmablas_ztranspose( AT(0, i), ldda, dAp, ldda, m,  sb);
            
            // make sure that gpu queue is empty
            //cuCtxSynchronize();
            
            // do the small non-parallel computations
            if ( s > (i+1) ) {
#ifndef WITHOUTTRTRI
                 cublasZtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                              sb, sb, 
                              c_one, L2(i),      lddl,
                                     UT(i, i+1), lddu);
#else
                 cublasZtrsm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                              sb, sb, 
                              c_one, L(i),      lddl,
                                     UT(i, i+1), lddu);
#endif
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             sb, m, sb,
                             c_neg_one, UT(i, i+1), lddu, 
                                        AT(0, i  ), ldda,
                             c_one,     AT(0, i+1), ldda );
            }
            else {
#ifndef WITHOUTTRTRI
                cublasZtrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n-mindim, sb, 
                             c_one, L2(i),      lddl,
                                    UT(i, i+1), lddu);
#else
                cublasZtrsm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n-mindim, sb, 
                             c_one, L(i),      lddl,
                                    UT(i, i+1), lddu);
#endif
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-mindim, m, sb,
                             c_neg_one, UT(i, i+1), lddu, 
                                        AT(0, i  ), ldda,
                             c_one,     AT(0, i+1), ldda );
            }
        }

	if ( (storev == 'C') || (storev == 'c') ) {
            magmablas_zgetmo_out( dU, dUT, lddu, m,  n );
            magmablas_zgetmo_out( dA, dAT, ldda, m,  n );
        }
    }
    return MAGMA_SUCCESS;
}
