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
magma_zunmqr_gpu(char side, char trans,
                 magma_int_t m, magma_int_t n, magma_int_t k,
		 cuDoubleComplex *dA,    magma_int_t ldda, 
                 cuDoubleComplex *tau,
                 cuDoubleComplex *dC,    magma_int_t lddc,
		 cuDoubleComplex *hwork, magma_int_t lwork,
                 cuDoubleComplex *dT,    magma_int_t nb, 
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

    A       (input) COMPLEX_16 array, dimension (LDDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            ZGEQRF in the first k columns of its array argument A.
            A is modified by the routine but restored on exit.

    LDDA     (input) INTEGER
            The leading dimension of the array A.
            If SIDE = 'L', LDDA >= max(1,M);
            if SIDE = 'R', LDDA >= max(1,N).

    TAU     (input) COMPLEX_16 array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZGEQRF.

    C       (input/output) COMPLEX_16 array, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q\*\*H*C or C*Q\*\*H or C*Q.

    LDDC     (input) INTEGER
            The leading dimension of the array C. LDDC >= max(1,M).

    HWORK    (workspace/output) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, HWORK(1) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The dimension of the array HWORK.
            If SIDE = 'L', LWORK >= max(1,N);
            if SIDE = 'R', LWORK >= max(1,M).
            For optimum performance LWORK >= N*NB if SIDE = 'L', and
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal
            blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the HWORK array, and no error
            message related to LWORK is issued by XERBLA.

    DT      (input) COMPLEX_16 array that is the output (the 9th argument)
            of magma_zgeqrf_gpu.

    NB      (input) INTEGER
            This is the blocking size that was used in pre-computing DT, e.g.,
            the blocking size used in magma_zgeqrf_gpu.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */

    #define a_ref(a_1,a_2) (dA+(a_2)*(ldda) + (a_1))
    #define c_ref(a_1,a_2) (dC+(a_2)*(lddc) + (a_1))
    #define t_ref(a_1)     (dT+(a_1)*nb)

    cuDoubleComplex c_one = MAGMA_Z_ONE;

    char side_[2] = {side, 0};
    char trans_[2] = {trans, 0};

    cuDoubleComplex *dwork;
    magma_int_t i, lddwork;

    magma_int_t i1, i2, i3, ib, ic, jc, mi, ni, nq, nw, ret;
    long int left, notran, lquery;
    static magma_int_t lwkopt;

    /* Function Body */
    *info = 0;
    left   = lapackf77_lsame(side_, "L");
    notran = lapackf77_lsame(trans_, "N");
    lquery = (lwork == -1);

    if (!left || notran)
      printf("zunmqr_gpu called with arguments not yet supported\n");

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
	nq = m;
	nw = n;
    } else {
	nq = n;
	nw = m;
    }
    if ( (!left) && (!lapackf77_lsame(side_, "R")) ) {
	*info = -1;
    } else if ( (!notran) && (!lapackf77_lsame(trans_, MagmaConjTransStr)) ) {
	*info = -2;
    } else if (m < 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (k < 0 || k > nq) {
	*info = -5;
    } else if (ldda < max(1,nq)) {
	*info = -7;
    } else if (lddc < max(1,m)) {
	*info = -10;
    } else if (lwork < max(1,nw) && ! lquery) {
	*info = -12;
    }

    lwkopt = (abs(m-k) + nb + 2*(n))*nb;
    hwork[0] = MAGMA_Z_MAKE( lwkopt, 0 );

    if (*info != 0) {
	return MAGMA_ERR_ILLEGAL_VALUE;
    } else if (lquery) {
	return MAGMA_SUCCESS;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
	hwork[0] = c_one;
	return MAGMA_SUCCESS;
    }

    lddwork= k;
    dwork  = dT+2*lddwork*nb;

    if ( (left && (! notran)) || ( (!left) && notran ) ) {
        i1 = 0;
        i2 = k-nb;
        i3 = nb;
    } else {
        i1 = (k - 1 - nb) / nb * nb;
        i2 = 0;
        i3 = -nb;
    }

    if (left) {
        ni = n;
        jc = 0;
    } else {
        mi = m;
        ic = 0;
    }

    if (nb < k)
    {
        for (i=i1; i3<0 ? i>i2 : i<i2; i+=i3)
        {
            ib = min(nb, k - i);
            if (left){
                mi = m - i;
                ic = i;
            }
            else {
                ni = n - i;
                jc = i;
            }
            ret = magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
				    mi, ni, ib, 
				    a_ref(i,  i ), ldda, t_ref(i), nb, 
				    c_ref(ic, jc), lddc, dwork,    nw);
	    if ( ret != MAGMA_SUCCESS )
	      return ret;
        }
    }
    else
    {
        i = i1;
    }

    /* Use unblocked code to multiply the last or only block. */
    if (i < k) {
        ib   = k-i;
        if (left){
            mi = m - i;
            ic = i;
        }
        else {
            ni = n - i;
            jc = i;
        }

        cublasGetMatrix(mi, ib, sizeof(cuDoubleComplex), 
			a_ref(i,  i ), ldda, hwork, mi);
        cublasGetMatrix(mi, ni, sizeof(cuDoubleComplex), 
			c_ref(ic, jc), lddc, hwork+mi*ib, mi);

        magma_int_t lhwork = lwork - mi*(ib + ni);
        lapackf77_zunmqr( MagmaLeftStr, MagmaConjTransStr, 
                          &mi, &ni, &ib, 
                          hwork,       &mi, tau+i, 
                          hwork+mi*ib, &mi, 
                          hwork+mi*(ib+ni), &lhwork, info);

        // send the updated part of c back to the GPU
        cublasSetMatrix(mi, ni, sizeof(cuDoubleComplex),
                        hwork+mi*ib, mi, c_ref(ic, jc), lddc);
    }

    return MAGMA_SUCCESS;
    /* End of MAGMA_ZUNMQR_GPU */
}
