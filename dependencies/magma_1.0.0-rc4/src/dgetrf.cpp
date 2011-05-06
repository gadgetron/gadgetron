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
  #define cublasDgemm magmablas_dgemm
  #define cublasDtrsm magmablas_dtrsm
#endif
// === End defining what BLAS to use =======================================


extern "C" magma_int_t
magma_dgetrf(magma_int_t m, magma_int_t n, double *a, magma_int_t lda, 
	     magma_int_t *ipiv, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    DGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.  This version does not
    require work space on the GPU passed as input. GPU memory is allocated
    in the routine.

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

    A       (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using cudaMallocHost.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -7, the GPU memory allocation failed
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================    */

#define inAT(i,j) (dAT + (i)*nb*ldda + (j)*nb)

    double *dAT, *dA, *da, *work;
    double c_one     = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;
    magma_int_t     iinfo, nb;

    /* Function Body */
    *info = 0;

    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,m))
        *info = -4;

    if (*info != 0)
        return MAGMA_ERR_ILLEGAL_VALUE;

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return MAGMA_SUCCESS;

    nb = magma_get_dgetrf_nb(m);

    if ( (nb <= 1) || (nb >= min(m,n)) ) {
        /* Use CPU code. */
	lapackf77_dgetrf(&m, &n, a, &lda, ipiv, info);
    } else {
        /* Use hybrid blocked code. */
        magma_int_t maxm, maxn, ldda, maxdim;
        magma_int_t i, rows, cols, s = min(m, n)/nb;

        maxm = ((m + 31)/32)*32;
        maxn = ((n + 31)/32)*32;
        maxdim = max(maxm, maxn);

        ldda = maxn;
        work = a;

        if (maxdim*maxdim < 2*maxm*maxn)
        {
	    if (CUBLAS_STATUS_SUCCESS != cublasAlloc(nb*maxm+maxdim*maxdim, sizeof(double), (void**)&dA) ) {
	        *info = -7;
		return MAGMA_ERR_CUBLASALLOC;
	    }
	    da = dA + nb*maxm;
	    
	    ldda = maxdim;
	    cublasSetMatrix( m, n, sizeof(double), a, lda, da, ldda);
	    
	    dAT = da;
	    magmablas_dinplace_transpose( dAT, ldda, ldda );
	}
        else
        {
	    if (CUBLAS_STATUS_SUCCESS != cublasAlloc((nb+maxn)*maxm, sizeof(double), (void**)&dA) ) {
	        *info = -7;
		return MAGMA_ERR_CUBLASALLOC;
	    }
	    da = dA + nb*maxm;
	    
	    cublasSetMatrix( m, n, sizeof(double), a, lda, da, maxm);
	    
	    if (CUBLAS_STATUS_SUCCESS != cublasAlloc(maxm*maxn, sizeof(double), (void**)&dAT) ) {
		cublasFree(dA);
		return MAGMA_ERR_CUBLASALLOC;
	    }

	    magmablas_dtranspose2( dAT, ldda, da, maxm, m, n );
	}
	
        lapackf77_dgetrf( &m, &nb, work, &lda, ipiv, &iinfo);

        for( i = 0; i < s; i++ )
        {
	    // download i-th panel
	    cols = maxm - i*nb;
	    
	    if (i>0){
		magmablas_dtranspose( dA, cols, inAT(i,i), ldda, nb, cols );
		cublasGetMatrix( m-i*nb, nb, sizeof(double), dA, cols, work, lda);
		
		// make sure that gpu queue is empty
		cuCtxSynchronize();
		
		cublasDtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
			     n - (i+1)*nb, nb, 
			     c_one, inAT(i-1,i-1), ldda, 
			            inAT(i-1,i+1), ldda );
		cublasDgemm( MagmaNoTrans, MagmaNoTrans, 
			     n-(i+1)*nb, m-i*nb, nb, 
			     c_neg_one, inAT(i-1,i+1), ldda, 
			                inAT(i,  i-1), ldda, 
			     c_one,     inAT(i,  i+1), ldda );

		// do the cpu part
		rows = m - i*nb;
		lapackf77_dgetrf( &rows, &nb, work, &lda, ipiv+i*nb, &iinfo);
	    }
	    if (*info == 0 && iinfo > 0)
		*info = iinfo + i*nb;
	    magmablas_dpermute_long2( dAT, ldda, ipiv, nb, i*nb );

	    // upload i-th panel
	    cublasSetMatrix( m-i*nb, nb, sizeof(double), work, lda, dA, cols);
	    magmablas_dtranspose( inAT(i,i), ldda, dA, cols, cols, nb);

	    // do the small non-parallel computations
	    if (s > (i+1)){
		cublasDtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
			     nb, nb, 
			     c_one, inAT(i, i  ), ldda,
                                    inAT(i, i+1), ldda);
		cublasDgemm( MagmaNoTrans, MagmaNoTrans, 
			     nb, m-(i+1)*nb, nb, 
			     c_neg_one, inAT(i,   i+1), ldda,
                                        inAT(i+1, i  ), ldda, 
			     c_one,     inAT(i+1, i+1), ldda );
	    }
	    else{
		cublasDtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n-s*nb, nb,
                             c_one, inAT(i, i  ), ldda,
                                    inAT(i, i+1), ldda);
		cublasDgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-(i+1)*nb, m-(i+1)*nb, nb,
                             c_neg_one, inAT(i,   i+1), ldda,
                                        inAT(i+1, i  ), ldda, 
                             c_one,     inAT(i+1, i+1), ldda );
	    }
	}
	
        magma_int_t nb0 = min(m - s*nb, n - s*nb);
        rows = m - s*nb;
        cols = maxm - s*nb;

        magmablas_dtranspose2( dA, cols, inAT(s,s), ldda, nb0, rows);
        cublasGetMatrix(rows, nb0, sizeof(double), dA, cols, work, lda);

        // make sure that gpu queue is empty
        cuCtxSynchronize();

        // do the cpu part
        lapackf77_dgetrf( &rows, &nb0, work, &lda, ipiv+s*nb, &iinfo);
        if (*info == 0 && iinfo > 0)
            *info = iinfo + s*nb;
        magmablas_dpermute_long2( dAT, ldda, ipiv, nb0, s*nb );

        cublasSetMatrix(rows, nb0, sizeof(double), work, lda, dA, cols);
        magmablas_dtranspose2( inAT(s,s), ldda, dA, cols, rows, nb0);

        cublasDtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                     n-s*nb-nb0, nb0,
                     c_one, inAT(s, s),     ldda, 
                            inAT(s, s)+nb0, ldda);

        if (maxdim*maxdim< 2*maxm*maxn){
            magmablas_dinplace_transpose( dAT, ldda, ldda );
            cublasGetMatrix( m, n, sizeof(double), da, ldda, a, lda);
        } else {
            magmablas_dtranspose2( da, maxm, dAT, ldda, n, m );
            cublasGetMatrix( m, n, sizeof(double), da, maxm, a, lda);
            cublasFree(dAT);
        }

        cublasFree(dA);
    }
    
    return 0;
} /* magma_dgetrf */

#undef inAT
