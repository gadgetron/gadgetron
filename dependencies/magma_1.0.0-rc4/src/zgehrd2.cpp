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
magma_zgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi, 
	      cuDoubleComplex *a, magma_int_t lda,
	      cuDoubleComplex *tau, cuDoubleComplex *work, 
	      magma_int_t *lwork, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   
    DGEHRD reduces a COMPLEX_16 general matrix A to upper Hessenberg form H by   
    an orthogonal similarity transformation:  Q' * A * Q = H .   

    Arguments   
    =========   
    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows   
            and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally   
            set by a previous call to DGEBAL; otherwise they should be   
            set to 1 and N respectively. See Further Details.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)   
            On entry, the N-by-N general matrix to be reduced.   
            On exit, the upper triangle and the first subdiagonal of A   
            are overwritten with the upper Hessenberg matrix H, and the   
            elements below the first subdiagonal, with the array TAU,   
            represent the orthogonal matrix Q as a product of elementary   
            reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    TAU     (output) COMPLEX_16 array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further   
            Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to   
            zero.   

    WORK    (workspace/output) COMPLEX_16 array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   
    The matrix Q is represented as a product of (ihi-ilo) elementary   
    reflectors   

       Q = H(ilo) H(ilo+1) . . . H(ihi-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on   
    exit in A(i+2:ihi,i), and tau in TAU(i).   

    The contents of A are illustrated by the following example, with   
    n = 7, ilo = 2 and ihi = 6:   

    on entry,                        on exit,   

    ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )   
    (     a   a   a   a   a   a )    (      a   h   h   h   h   a )   
    (     a   a   a   a   a   a )    (      h   h   h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )   
    (                         a )    (                          a )   

    where a denotes an element of the original matrix A, h denotes a   
    modified element of the upper Hessenberg matrix H, and vi denotes an   
    element of the vector defining H(i).   

    This implementation follows the hybrid algorithm and notations described in

    S. Tomov and J. Dongarra, "Accelerating the reduction to upper Hessenberg
    form through hybrid GPU-based computing," University of Tennessee Computer
    Science Technical Report, UT-CS-09-642 (also LAPACK Working Note 219),
    May 24, 2009.
    =====================================================================    */


    cuDoubleComplex c_one = MAGMA_Z_ONE;
    cuDoubleComplex c_zero = MAGMA_Z_ZERO;

    magma_int_t nb = magma_get_zgehrd_nb(n);
    magma_int_t N = n, ldda = n;

    magma_int_t ib;
    magma_int_t nh, iws;
    magma_int_t nbmin, iinfo;
    magma_int_t ldwork;
    magma_int_t lquery;

    --tau;

    *info = 0;
    MAGMA_Z_SET2REAL( work[0], (double) n * nb );

    lquery = *lwork == -1;
    if (n < 0) {
	*info = -1;
    } else if (ilo < 1 || ilo > max(1,n)) {
	*info = -2;
    } else if (ihi < min(ilo,n) || ihi > n) {
	*info = -3;
    } else if (lda < max(1,n)) {
	*info = -5;
    } else if (*lwork < max(1,n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0)
      return 0;
    else if (lquery)
      return 0;

    /* Quick return if possible */
    nh = ihi - ilo + 1;
    if (nh <= 1) {
      work[0] = c_one;
      return 0;
    }

    cuDoubleComplex *da;
    cublasStatus status;
    status = cublasAlloc(N*ldda+2*N*nb+nb*nb, sizeof(cuDoubleComplex), (void**)&da);
    if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device memory allocation error (d_A)\n");
    }
    
    cuDoubleComplex *d_A    = da;
    cuDoubleComplex *d_work = da + (N+nb)*ldda;

    magma_int_t i__;

    cuDoubleComplex *t, *d_t;
    t   = (cuDoubleComplex *)malloc(nb*nb*sizeof(cuDoubleComplex));
    d_t = d_work + nb * ldda;

    zzero_nbxnb_block(nb, d_A+N*ldda, ldda);

    /* Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */
    for (i__ = 1; i__ < ilo; ++i__)
      tau[i__] = c_zero;
   
    for (i__ = max(1,ihi); i__ < n; ++i__)
      tau[i__] = c_zero;

    for(i__=0; i__< nb*nb; i__+=4)
      t[i__] = t[i__+1] = t[i__+2] = t[i__+3] = c_zero;

    nbmin = 2;
    iws = 1;
    if (nb > 1 && nb < nh) {

      /*  Determine when to cross over from blocked to unblocked code   
          (last block is always handled by unblocked code)              */
      if (nb < nh) {

	/* Determine if workspace is large enough for blocked code      */
	iws = n * nb;
	if (*lwork < iws) {

	  /*    Not enough workspace to use optimal NB:  determine the   
                minimum value of NB, and reduce NB or force use of   
                unblocked code                                          */
	  nbmin = nb;
	  if (*lwork >= n * nbmin)
	    nb = *lwork / n;
	  else 
	    nb = 1;
	}
      }
    }
    ldwork = n;

    if (nb < nbmin || nb >= nh) {
      /* Use unblocked code below */
      i__ = ilo;
    } else {

      /* Use blocked code */

      /* Copy the matrix to the GPU */
      cublasSetMatrix(N, N, sizeof(cuDoubleComplex), a+(ilo-1)*(lda),lda, d_A, ldda);

      for (i__ = ilo; i__ < ihi - nb; i__ += nb) {
	/* Computing MIN */
	ib = min(nb, ihi - i__);

	/*   Reduce columns i:i+ib-1 to Hessenberg form, returning the   
             matrices V and T of the block reflector H = I - V*T*V'   
             which performs the reduction, and also the matrix Y = A*V*T */

	/*   Get the current panel (no need for the 1st iteration) */
	cublasGetMatrix(ihi-i__+ilo, ib, sizeof(cuDoubleComplex), 
			d_A + (i__ - ilo)*ldda  + i__ - ilo, ldda,
			a   + (i__ -   1 )*(lda)+ i__ - 1   , lda);
	
	magma_zlahr2(ihi, i__, ib, 
		     d_A + (i__ - ilo)*ldda, 
		     d_A + N*ldda + 1,
		     a   + (i__ -   1 )*(lda) , lda, 
		     &tau[i__], t, nb, work, ldwork);

	/* Copy T from the CPU to D_T on the GPU */
	cublasSetMatrix(nb, nb, sizeof(cuDoubleComplex), t, nb, d_t, nb);

	magma_zlahru(ihi, i__ - ilo, ib, 
		     a   + (i__ -   1 )*(lda), lda,
		     d_A + (i__ - ilo)*ldda, 
		     d_A + (i__ - ilo)*ldda + i__ - 1,
		     d_A + N*ldda, d_t, d_work);
      }
    }

    /* Use unblocked code to reduce the rest of the matrix */
    if (!(nb < nbmin || nb >= nh))
      cublasGetMatrix(n, n-i__+1, sizeof(cuDoubleComplex), 
		      d_A+ (i__-1)*ldda, ldda, 
		      a  + (i__-1)*(lda), lda);
    lapackf77_zgehd2(&n, &i__, &ihi, a, &lda, &tau[1], work, &iinfo);
    MAGMA_Z_SET2REAL( work[0], (double) iws );
    
    cublasFree(da);
    free(t);
 
    return 0;
} /* magma_zgehrd2 */

