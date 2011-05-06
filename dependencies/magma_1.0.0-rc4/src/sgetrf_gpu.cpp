/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated s

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_s
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasSgemm magmablas_sgemm
  #define cublasStrsm magmablas_strsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_sgetrf_gpu(magma_int_t m, magma_int_t n, 
		 float *dA, magma_int_t ldda,
		 magma_int_t *ipiv, magma_int_t *info)
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

    A       (input/output) REAL array on the GPU, dimension (LDDA,N).
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

    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;

    magma_int_t iinfo, nb;
    magma_int_t maxm, maxn, mindim;
    magma_int_t i, rows, cols, s, lddat, lddwork;
    float *dAT, *dAP, *work;

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
    nb     = magma_get_sgetrf_nb(m);
    s      = mindim / nb;

    if (nb <= 1 || nb >= min(m,n)) {
	/* Use CPU code. */
	work = (float*)malloc(m * n * sizeof(float));
	cublasGetMatrix(m, n, sizeof(float), dA, ldda, work, m);
	lapackf77_sgetrf(&m, &n, work, &m, ipiv, info);
	cublasSetMatrix(m, n, sizeof(float), work, m, dA, ldda);
	free(work);
    }
    else {
	/* Use hybrid blocked code. */
	maxm = ((m + 31)/32)*32;
	maxn = ((n + 31)/32)*32;

	lddat   = maxn;
	lddwork = maxm;

	dAT = dA;

	if ( CUBLAS_STATUS_SUCCESS != cublasAlloc(nb*maxm, sizeof(float), (void**)&dAP) ) {
	    return MAGMA_ERR_CUBLASALLOC;
	}

	if ((m == n) && (m % 32 == 0) && (ldda%32 == 0))
	    magmablas_sinplace_transpose( dAT, ldda, lddat );
	else {
	    if ( CUBLAS_STATUS_SUCCESS != cublasAlloc(maxm*maxn, sizeof(float), (void**)&dAT) ) {
		cublasFree( dAP );
		return MAGMA_ERR_CUBLASALLOC;
	    }
	    magmablas_stranspose2( dAT, lddat, dA, ldda, m, n );
	}

	if ( cudaSuccess != cudaMallocHost( (void**)&work, maxm*nb*sizeof(float) ) ) {
	    cublasFree( dAP );
	    if (! ((m == n) && (m % 32 == 0) && (ldda%32 == 0)) )
		cublasFree( dAT );
	    return MAGMA_ERR_HOSTALLOC;
	}

	for( i=0; i<s; i++ )
            {
                // download i-th panel
                cols = maxm - i*nb;
                magmablas_stranspose( dAP, cols, inAT(i,i), lddat, nb, cols );
                cublasGetMatrix( m-i*nb, nb, sizeof(float), dAP, cols, work, lddwork);

                // make sure that gpu queue is empty
                cuCtxSynchronize();

                if ( i>0 ){
                    cublasStrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 n - (i+1)*nb, nb, 
                                 c_one, inAT(i-1,i-1), lddat, 
                                        inAT(i-1,i+1), lddat );
                    cublasSgemm( MagmaNoTrans, MagmaNoTrans, 
                                 n-(i+1)*nb, m-i*nb, nb, 
                                 c_neg_one, inAT(i-1,i+1), lddat, 
                                            inAT(i,  i-1), lddat, 
                                 c_one,     inAT(i,  i+1), lddat );
                }

                // do the cpu part
                rows = m - i*nb;
                lapackf77_sgetrf( &rows, &nb, work, &lddwork, ipiv+i*nb, &iinfo);
                if ( (*info == 0) && (iinfo > 0) )
                    *info = iinfo + i*nb;

                magmablas_spermute_long2( dAT, lddat, ipiv, nb, i*nb );

                // upload i-th panel
                cublasSetMatrix(m-i*nb, nb, sizeof(float), work, lddwork, dAP, maxm);
                magmablas_stranspose(inAT(i,i), lddat, dAP, maxm, cols, nb);

                // do the small non-parallel computations
                if ( s > (i+1) ) {
                    cublasStrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 nb, nb, 
                                 c_one, inAT(i, i  ), lddat,
                                        inAT(i, i+1), lddat);
                    cublasSgemm( MagmaNoTrans, MagmaNoTrans, 
                                 nb, m-(i+1)*nb, nb, 
                                 c_neg_one, inAT(i,   i+1), lddat,
                                            inAT(i+1, i  ), lddat, 
                                 c_one,     inAT(i+1, i+1), lddat );
                }
                else {
                    cublasStrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 n-s*nb, nb, 
                                 c_one, inAT(i, i  ), lddat,
                                        inAT(i, i+1), lddat);
                    cublasSgemm( MagmaNoTrans, MagmaNoTrans, 
                                 n-(i+1)*nb, m-(i+1)*nb, nb,
                                 c_neg_one, inAT(i,   i+1), lddat,
                                            inAT(i+1, i  ), lddat, 
                                 c_one,     inAT(i+1, i+1), lddat );
                }
            }

	magma_int_t nb0 = min(m - s*nb, n - s*nb);
	rows = m - s*nb;
	cols = maxm - s*nb;

	magmablas_stranspose2( dAP, maxm, inAT(s,s), lddat, nb0, rows);
	cublasGetMatrix(rows, nb0, sizeof(float), dAP, maxm, work, lddwork);

	// make sure that gpu queue is empty
	cuCtxSynchronize();

	// do the cpu part
	lapackf77_sgetrf( &rows, &nb0, work, &lddwork, ipiv+s*nb, &iinfo);
	if ( (*info == 0) && (iinfo > 0) )
	    *info = iinfo + s*nb;
	magmablas_spermute_long2( dAT, lddat, ipiv, nb0, s*nb );

	// upload i-th panel
	cublasSetMatrix(rows, nb0, sizeof(float), work, lddwork, dAP, maxm);
	magmablas_stranspose2( inAT(s,s), lddat, dAP, maxm, rows, nb0);

	cublasStrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                     n-s*nb-nb0, nb0,
		     c_one, inAT(s,s),     lddat, 
                            inAT(s,s)+nb0, lddat);

	if ((m == n) && (m % 32 == 0) && (ldda%32 == 0))
	    magmablas_sinplace_transpose( dAT, lddat, ldda );
	else {
	    magmablas_stranspose2( dA, ldda, dAT, lddat, n, m );
	    cublasFree(dAT);
	}

	cublasFree(dAP);
	cudaFreeHost(work);
    }
    return MAGMA_SUCCESS;

    /* End of MAGMA_SGETRF_GPU */
}

#undef inAT
