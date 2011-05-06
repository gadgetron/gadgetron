/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/
#include "common_magma.h"

extern "C" int
magma_zunmtr(char side, char uplo, char trans,
             int m, int n, 
	     cuDoubleComplex *a,    int lda, 
             cuDoubleComplex *tau, 
             cuDoubleComplex *c,    int ldc,
	     cuDoubleComplex *work, int lwork, 
             int *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   

    ZUNMTR overwrites the general complex M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q\*\*H * C       C * Q\*\*H   

    where Q is a complex orthogonal matrix of order nq, with nq = m if   
    SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of   
    nq-1 elementary reflectors, as returned by SSYTRD:   

    if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);   

    if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q\*\*H from the Left;   
            = 'R': apply Q or Q\*\*H from the Right.   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangle of A contains elementary reflectors   
                   from SSYTRD;   
            = 'L': Lower triangle of A contains elementary reflectors   
                   from SSYTRD.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q\*\*H.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    A       (input) COMPLEX_16 array, dimension   
                                 (LDA,M) if SIDE = 'L'   
                                 (LDA,N) if SIDE = 'R'   
            The vectors which define the elementary reflectors, as   
            returned by SSYTRD.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.   

    TAU     (input) COMPLEX_16 array, dimension   
                                 (M-1) if SIDE = 'L'   
                                 (N-1) if SIDE = 'R'   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by SSYTRD.   

    C       (input/output) COMPLEX_16 array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by Q*C or Q\*\*H*C or C*Q\*\*H or C*Q.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace/output) COMPLEX_16 array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If SIDE = 'L', LWORK >= max(1,N);   
            if SIDE = 'R', LWORK >= max(1,M).   
            For optimum performance LWORK >= N*NB if SIDE = 'L', and   
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal   
            blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================    */

   
    cuDoubleComplex c_one = MAGMA_Z_ONE;

    char side_[2]  = {side, 0};
    char uplo_[2]  = {uplo, 0};
    char trans_[2] = {trans, 0};
    int a_dim1, a_offset, c_dim1, c_offset, i__2;
    static int i1, i2, nb, mi, ni, nq, nw;
    long int left, upper, lquery;
    static int iinfo;
    static int lwkopt;

    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = ldc;
    c_offset = 1 + c_dim1;
    c -= c_offset;

    /* Function Body */
    *info = 0;
    left   = lapackf77_lsame(side_, "L");
    upper  = lapackf77_lsame(uplo_, "U");
    lquery = lwork == -1;

    /* NQ is the order of Q and NW is the minimum dimension of WORK */

    if (left) {
	nq = m;
	nw = n;
    } else {
	nq = n;
	nw = m;
    }
    if (! left && ! lapackf77_lsame(side_, "R")) {
	*info = -1;
    } else if (! upper && ! lapackf77_lsame(uplo_, "L")) {
	*info = -2;
    } else if (! lapackf77_lsame(trans_, "N") && 
               ! lapackf77_lsame(trans_, "C")) {
	*info = -3;
    } else if (m < 0) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max(1,nq)) {
	*info = -7;
    } else if (ldc < max(1,m)) {
	*info = -10;
    } else if (lwork < max(1,nw) && ! lquery) {
	*info = -12;
    }

    if (*info == 0) 
      {
	nb = 32;
	lwkopt = max(1,nw) * nb;
	MAGMA_Z_SET2REAL( work[0], lwkopt );
      }

    if (*info != 0) {
	i__2 = -(*info);
	return MAGMA_ERR_ILLEGAL_VALUE;
    } else if (lquery) {
	return MAGMA_SUCCESS;
    }

/*     Quick return if possible */

    if (m == 0 || n == 0 || nq == 1) {
	work[0] = c_one;
	return MAGMA_SUCCESS;
    }

    if (left) {
	mi = m - 1;
	ni = n;
    } else {
	mi = m;
	ni = n - 1;
    }

    if (upper) 
      {
	/* Q was determined by a call to SSYTRD with UPLO = 'U' */
	i__2 = nq - 1;
	lapackf77_zunmql(side_, trans_, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], &lda, 
                         tau, &c[c_offset], &ldc, work, &lwork, &iinfo);
      }
    else 
      {
	/* Q was determined by a call to SSYTRD with UPLO = 'L' */
	if (left) {
	    i1 = 2;
	    i2 = 1;
	} else {
	    i1 = 1;
	    i2 = 2;
	}
	i__2 = nq - 1;
	magma_zunmqr(side, trans, mi, ni, i__2, &a[a_dim1 + 2], lda, tau,
		     &c[i1 + i2 * c_dim1], ldc, work, lwork, &iinfo);
      }
    MAGMA_Z_SET2REAL( work[0], lwkopt );
    return MAGMA_SUCCESS;
/*     End of ZUNMTR */

} /* zunmtr_ */

