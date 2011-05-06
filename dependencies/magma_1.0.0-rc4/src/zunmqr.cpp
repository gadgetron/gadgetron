/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_zunmqr(const char side, const char trans, 
             magma_int_t m, magma_int_t n, magma_int_t k, 
             cuDoubleComplex *a,    magma_int_t lda, 
             cuDoubleComplex *tau, 
             cuDoubleComplex *c,    magma_int_t ldc,
	     cuDoubleComplex *work, magma_int_t lwork, 
             magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   
    ZUNMQR overwrites the general complex M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q\*\*H * C       C * Q\*\*H   

    where Q is a complex orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(1) H(2) . . . H(k)   

    as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q\*\*H from the Left;   
            = 'R': apply Q or Q\*\*H from the Right.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q\*\*H.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) COMPLEX_16 array, dimension (LDA,K)   
            The i-th column must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by   
            ZGEQRF in the first k columns of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            If SIDE = 'L', LDA >= max(1,M);   
            if SIDE = 'R', LDA >= max(1,N).   

    TAU     (input) COMPLEX_16 array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by ZGEQRF.   

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
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   */

    
    cuDoubleComplex c_one = MAGMA_Z_ONE;

    char side_[2] = {side, 0};
    char trans_[2] = {trans, 0};

    // TTT --------------------------------------------------------------------
    cuDoubleComplex *dwork, *dc;
    cublasAlloc((m)*(n), sizeof(cuDoubleComplex), (void**)&dc);
    cublasAlloc(2*(m+64)*64, sizeof(cuDoubleComplex), (void**)&dwork);
    
    cublasSetMatrix( m, n, sizeof(cuDoubleComplex), c, ldc, dc, ldc);
    dc -= (1 + m);
    //-------------------------------------------------------------------------

    int a_dim1, a_offset, c_dim1, c_offset, i__4, i__5;
    /* Local variables */
    static int i__;
    static cuDoubleComplex t[2*4160]	/* was [65][64] */;
    static int i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    long int left, notran, lquery;
    static int nbmin, iinfo;
    static int ldwork, lwkopt;

    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = ldc;
    c_offset = 1 + c_dim1;
    c -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lapackf77_lsame(side_, "L");
    notran = lapackf77_lsame(trans_, "N");
    lquery = (lwork == -1);

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
    } else if (! notran && ! lapackf77_lsame(trans_, "T")) {
	*info = -2;
    } else if (m < 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (k < 0 || k > nq) {
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
	/* Determine the block size.  NB may be at most NBMAX, where NBMAX   
	   is used to define the local array T.    */
	nb = 64;
	lwkopt = max(1,nw) * nb;
	MAGMA_Z_SET2REAL( work[1], lwkopt );
    }

    if (*info != 0) {
	return 0;
    } else if (lquery) {
	return 0;
    }

    /* Quick return if possible */

    if (m == 0 || n == 0 || k == 0) {
	work[1] = c_one;
	return 0;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < k) {
	iws = nw * nb;
	if (lwork < iws) {
	    nb = lwork / ldwork;
	    nbmin = 64;
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= k) 
      {
	/* Use unblocked code */
	lapackf77_zunm2r(side_, trans_, &m, &n, &k, &a[a_offset], &lda, &tau[1], 
		&c[c_offset], &ldc, &work[1], &iinfo);
      } 
    else 
      {
	/* Use blocked code */
	if ( ( left && (! notran) ) ||  ( (! left) && notran ) ) {
	    i1 = 1;
	    i2 = k;
	    i3 = nb;
	} else {
	    i1 = (k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = n;
	    jc = 1;
	} else {
	    mi = m;
	    ic = 1;
	}
	
	for (i__ = i1; i3 < 0 ? i__ >= i2 : i__ <= i2; i__ += i3) 
	  {
	    /* Computing MIN */
	    i__4 = nb, i__5 = k - i__ + 1;
	    ib = min(i__4,i__5);

	    /* Form the triangular factor of the block reflector   
	       H = H(i) H(i+1) . . . H(i+ib-1) */
	    i__4 = nq - i__ + 1;
	    lapackf77_zlarft("F", "C", &i__4, &ib, &a[i__ + i__ * a_dim1], &lda, 
		    &tau[i__], t, &ib);

	    // TTT ------------------------------------------------------------
	    zpanel_to_q('U', ib, &a[i__ + i__ * a_dim1], lda, t+ib*ib);
	    cublasSetMatrix(i__4, ib, sizeof(cuDoubleComplex),
			    &a[i__ + i__ * a_dim1], lda, 
			    dwork, i__4);
	    zq_to_panel('U', ib, &a[i__ + i__ * a_dim1], lda, t+ib*ib);
	    //-----------------------------------------------------------------

	    if (left) 
	      {
		/* H or H' is applied to C(i:m,1:n) */
		mi = m - i__ + 1;
		ic = i__;
	      } 
	    else 
	      {
		/* H or H' is applied to C(1:m,i:n) */
		ni = n - i__ + 1;
		jc = i__;
	      }
	    
	    /* Apply H or H' */
	    // TTT ------------------------------------------------------------
	    //printf("%5d %5d %5d\n", mi, ni, ic + 1 + m);
	    cublasSetMatrix(ib, ib, sizeof(cuDoubleComplex), t, ib, dwork+i__4*ib, ib);
	    magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
			      mi, ni, ib,
			      dwork, i__4, dwork+i__4*ib, ib,
			      &dc[ic + jc * c_dim1], ldc, 
			      dwork+i__4*ib + ib*ib, ni);
	    //-----------------------------------------------------------------
	    /*
	    lapackf77_zlarfb(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, 
		    &a[i__ + i__ * a_dim1], &lda, t, &c__65, 
		    &c[ic + jc * c_dim1], &ldc, &work[1], &ldwork);
	    */
	  }
      }
    MAGMA_Z_SET2REAL( work[1], lwkopt );

    dc += (1 + m);
    cublasFree(dc);
    cublasFree(dwork);

    return 0;
    
/*     End of ZUNMQR */

} /* zunmqr_ */


