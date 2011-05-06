/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasZgemm magmablas_zgemm
  #define cublasZtrsm magmablas_ztrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_zgetrf_gpu(magma_int_t m, magma_int_t n, 
		 cuDoubleComplex *dA, magma_int_t ldda,
		 magma_int_t *ipiv, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    ZGETRF computes an LU factorization of a general M-by-N matrix A
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

    A       (input/output) COMPLEX_16 array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDDA     (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

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

#define inAT(i,j) (dAT + (i)*nb*lddat + (j)*nb)

    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    magma_int_t iinfo, nb;
    magma_int_t maxm, maxn, mindim;
    magma_int_t i, rows, cols, s, lddat, lddwork;
    cuDoubleComplex *dAT, *dAP, *work;

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
    nb     = magma_get_zgetrf_nb(m);
    s      = mindim / nb;

    if (nb <= 1 || nb >= min(m,n)) {
	/* Use CPU code. */
	work = (cuDoubleComplex*)malloc(m * n * sizeof(cuDoubleComplex));
	cublasGetMatrix(m, n, sizeof(cuDoubleComplex), dA, ldda, work, m);
	lapackf77_zgetrf(&m, &n, work, &m, ipiv, info);
	cublasSetMatrix(m, n, sizeof(cuDoubleComplex), work, m, dA, ldda);
	free(work);
    }
    else {
	/* Use hybrid blocked code. */
	maxm = ((m + 31)/32)*32;
	maxn = ((n + 31)/32)*32;

	lddat   = maxn;
	lddwork = maxm;

	dAT = dA;

	if ( CUBLAS_STATUS_SUCCESS != cublasAlloc(nb*maxm, sizeof(cuDoubleComplex), (void**)&dAP) ) {
	    return MAGMA_ERR_CUBLASALLOC;
	}

	if ((m == n) && (m % 32 == 0) && (ldda%32 == 0))
	    magmablas_zinplace_transpose( dAT, ldda, lddat );
	else {
	    if ( CUBLAS_STATUS_SUCCESS != cublasAlloc(maxm*maxn, sizeof(cuDoubleComplex), (void**)&dAT) ) {
		cublasFree( dAP );
		return MAGMA_ERR_CUBLASALLOC;
	    }
	    magmablas_ztranspose2( dAT, lddat, dA, ldda, m, n );
	}

	if ( cudaSuccess != cudaMallocHost( (void**)&work, maxm*nb*sizeof(cuDoubleComplex) ) ) {
	    cublasFree( dAP );
	    if (! ((m == n) && (m % 32 == 0) && (ldda%32 == 0)) )
		cublasFree( dAT );
	    return MAGMA_ERR_HOSTALLOC;
	}

	for( i=0; i<s; i++ )
            {
                // download i-th panel
                cols = maxm - i*nb;
                magmablas_ztranspose( dAP, cols, inAT(i,i), lddat, nb, cols );
                cublasGetMatrix( m-i*nb, nb, sizeof(cuDoubleComplex), dAP, cols, work, lddwork);

                // make sure that gpu queue is empty
                cuCtxSynchronize();

                if ( i>0 ){
                    cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 n - (i+1)*nb, nb, 
                                 c_one, inAT(i-1,i-1), lddat, 
                                        inAT(i-1,i+1), lddat );
                    cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                                 n-(i+1)*nb, m-i*nb, nb, 
                                 c_neg_one, inAT(i-1,i+1), lddat, 
                                            inAT(i,  i-1), lddat, 
                                 c_one,     inAT(i,  i+1), lddat );
                }

                // do the cpu part
                rows = m - i*nb;
                lapackf77_zgetrf( &rows, &nb, work, &lddwork, ipiv+i*nb, &iinfo);
                if ( (*info == 0) && (iinfo > 0) )
                    *info = iinfo + i*nb;

                magmablas_zpermute_long2( dAT, lddat, ipiv, nb, i*nb );

                // upload i-th panel
                cublasSetMatrix(m-i*nb, nb, sizeof(cuDoubleComplex), work, lddwork, dAP, maxm);
                magmablas_ztranspose(inAT(i,i), lddat, dAP, maxm, cols, nb);

                // do the small non-parallel computations
                if ( s > (i+1) ) {
                    cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 nb, nb, 
                                 c_one, inAT(i, i  ), lddat,
                                        inAT(i, i+1), lddat);
                    cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                                 nb, m-(i+1)*nb, nb, 
                                 c_neg_one, inAT(i,   i+1), lddat,
                                            inAT(i+1, i  ), lddat, 
                                 c_one,     inAT(i+1, i+1), lddat );
                }
                else {
                    cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 n-s*nb, nb, 
                                 c_one, inAT(i, i  ), lddat,
                                        inAT(i, i+1), lddat);
                    cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                                 n-(i+1)*nb, m-(i+1)*nb, nb,
                                 c_neg_one, inAT(i,   i+1), lddat,
                                            inAT(i+1, i  ), lddat, 
                                 c_one,     inAT(i+1, i+1), lddat );
                }
            }

	magma_int_t nb0 = min(m - s*nb, n - s*nb);
	rows = m - s*nb;
	cols = maxm - s*nb;

	magmablas_ztranspose2( dAP, maxm, inAT(s,s), lddat, nb0, rows);
	cublasGetMatrix(rows, nb0, sizeof(cuDoubleComplex), dAP, maxm, work, lddwork);

	// make sure that gpu queue is empty
	cuCtxSynchronize();

	// do the cpu part
	lapackf77_zgetrf( &rows, &nb0, work, &lddwork, ipiv+s*nb, &iinfo);
	if ( (*info == 0) && (iinfo > 0) )
	    *info = iinfo + s*nb;
	magmablas_zpermute_long2( dAT, lddat, ipiv, nb0, s*nb );

	// upload i-th panel
	cublasSetMatrix(rows, nb0, sizeof(cuDoubleComplex), work, lddwork, dAP, maxm);
	magmablas_ztranspose2( inAT(s,s), lddat, dAP, maxm, rows, nb0);

	cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                     n-s*nb-nb0, nb0,
		     c_one, inAT(s,s),     lddat, 
                            inAT(s,s)+nb0, lddat);

	if ((m == n) && (m % 32 == 0) && (ldda%32 == 0))
	    magmablas_zinplace_transpose( dAT, lddat, ldda );
	else {
	    magmablas_ztranspose2( dA, ldda, dAT, lddat, n, m );
	    cublasFree(dAT);
	}

	cublasFree(dAP);
	cudaFreeHost(work);
    }
    return MAGMA_SUCCESS;

    /* End of MAGMA_ZGETRF_GPU */
}

#undef inAT
