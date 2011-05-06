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
magma_zgetrf(magma_int_t m, magma_int_t n, cuDoubleComplex *a, magma_int_t lda, 
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

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)
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

    cuDoubleComplex *dAT, *dA, *da, *work;
    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
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

    nb = magma_get_zgetrf_nb(m);

    if ( (nb <= 1) || (nb >= min(m,n)) ) {
        /* Use CPU code. */
	lapackf77_zgetrf(&m, &n, a, &lda, ipiv, info);
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
	    if (CUBLAS_STATUS_SUCCESS != cublasAlloc(nb*maxm+maxdim*maxdim, sizeof(cuDoubleComplex), (void**)&dA) ) {
	        *info = -7;
		return MAGMA_ERR_CUBLASALLOC;
	    }
	    da = dA + nb*maxm;
	    
	    ldda = maxdim;
	    cublasSetMatrix( m, n, sizeof(cuDoubleComplex), a, lda, da, ldda);
	    
	    dAT = da;
	    magmablas_zinplace_transpose( dAT, ldda, ldda );
	}
        else
        {
	    if (CUBLAS_STATUS_SUCCESS != cublasAlloc((nb+maxn)*maxm, sizeof(cuDoubleComplex), (void**)&dA) ) {
	        *info = -7;
		return MAGMA_ERR_CUBLASALLOC;
	    }
	    da = dA + nb*maxm;
	    
	    cublasSetMatrix( m, n, sizeof(cuDoubleComplex), a, lda, da, maxm);
	    
	    if (CUBLAS_STATUS_SUCCESS != cublasAlloc(maxm*maxn, sizeof(cuDoubleComplex), (void**)&dAT) ) {
		cublasFree(dA);
		return MAGMA_ERR_CUBLASALLOC;
	    }

	    magmablas_ztranspose2( dAT, ldda, da, maxm, m, n );
	}
	
        lapackf77_zgetrf( &m, &nb, work, &lda, ipiv, &iinfo);

        for( i = 0; i < s; i++ )
        {
	    // download i-th panel
	    cols = maxm - i*nb;
	    
	    if (i>0){
		magmablas_ztranspose( dA, cols, inAT(i,i), ldda, nb, cols );
		cublasGetMatrix( m-i*nb, nb, sizeof(cuDoubleComplex), dA, cols, work, lda);
		
		// make sure that gpu queue is empty
		cuCtxSynchronize();
		
		cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
			     n - (i+1)*nb, nb, 
			     c_one, inAT(i-1,i-1), ldda, 
			            inAT(i-1,i+1), ldda );
		cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
			     n-(i+1)*nb, m-i*nb, nb, 
			     c_neg_one, inAT(i-1,i+1), ldda, 
			                inAT(i,  i-1), ldda, 
			     c_one,     inAT(i,  i+1), ldda );

		// do the cpu part
		rows = m - i*nb;
		lapackf77_zgetrf( &rows, &nb, work, &lda, ipiv+i*nb, &iinfo);
	    }
	    if (*info == 0 && iinfo > 0)
		*info = iinfo + i*nb;
	    magmablas_zpermute_long2( dAT, ldda, ipiv, nb, i*nb );

	    // upload i-th panel
	    cublasSetMatrix( m-i*nb, nb, sizeof(cuDoubleComplex), work, lda, dA, cols);
	    magmablas_ztranspose( inAT(i,i), ldda, dA, cols, cols, nb);

	    // do the small non-parallel computations
	    if (s > (i+1)){
		cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
			     nb, nb, 
			     c_one, inAT(i, i  ), ldda,
                                    inAT(i, i+1), ldda);
		cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
			     nb, m-(i+1)*nb, nb, 
			     c_neg_one, inAT(i,   i+1), ldda,
                                        inAT(i+1, i  ), ldda, 
			     c_one,     inAT(i+1, i+1), ldda );
	    }
	    else{
		cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n-s*nb, nb,
                             c_one, inAT(i, i  ), ldda,
                                    inAT(i, i+1), ldda);
		cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-(i+1)*nb, m-(i+1)*nb, nb,
                             c_neg_one, inAT(i,   i+1), ldda,
                                        inAT(i+1, i  ), ldda, 
                             c_one,     inAT(i+1, i+1), ldda );
	    }
	}
	
        magma_int_t nb0 = min(m - s*nb, n - s*nb);
        rows = m - s*nb;
        cols = maxm - s*nb;

        magmablas_ztranspose2( dA, cols, inAT(s,s), ldda, nb0, rows);
        cublasGetMatrix(rows, nb0, sizeof(cuDoubleComplex), dA, cols, work, lda);

        // make sure that gpu queue is empty
        cuCtxSynchronize();

        // do the cpu part
        lapackf77_zgetrf( &rows, &nb0, work, &lda, ipiv+s*nb, &iinfo);
        if (*info == 0 && iinfo > 0)
            *info = iinfo + s*nb;
        magmablas_zpermute_long2( dAT, ldda, ipiv, nb0, s*nb );

        cublasSetMatrix(rows, nb0, sizeof(cuDoubleComplex), work, lda, dA, cols);
        magmablas_ztranspose2( inAT(s,s), ldda, dA, cols, rows, nb0);

        cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                     n-s*nb-nb0, nb0,
                     c_one, inAT(s, s),     ldda, 
                            inAT(s, s)+nb0, ldda);

        if (maxdim*maxdim< 2*maxm*maxn){
            magmablas_zinplace_transpose( dAT, ldda, ldda );
            cublasGetMatrix( m, n, sizeof(cuDoubleComplex), da, ldda, a, lda);
        } else {
            magmablas_ztranspose2( da, maxm, dAT, ldda, n, m );
            cublasGetMatrix( m, n, sizeof(cuDoubleComplex), da, maxm, a, lda);
            cublasFree(dAT);
        }

        cublasFree(dA);
    }
    
    return 0;
} /* magma_zgetrf */

#undef inAT
