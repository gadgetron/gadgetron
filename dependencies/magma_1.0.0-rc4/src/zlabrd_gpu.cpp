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
//  #define cublasZgemv magmablas_zgemv
#endif
// === End defining what BLAS to use ======================================

extern "C" magma_int_t 
magma_zlabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb,
                  cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex *da, magma_int_t ldda,
                  double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,
                  cuDoubleComplex *x, magma_int_t ldx, cuDoubleComplex *dx, magma_int_t lddx,
                  cuDoubleComplex *y, magma_int_t ldy, cuDoubleComplex *dy, magma_int_t lddy)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   
    ZLABRD reduces the first NB rows and columns of a complex general   
    m by n matrix A to upper or lower bidiagonal form by an orthogonal   
    transformation Q' * A * P, and returns the matrices X and Y which   
    are needed to apply the transformation to the unreduced part of A.   

    If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower   
    bidiagonal form.   

    This is an auxiliary routine called by SGEBRD   

    Arguments   
    =========   
    M       (input) INTEGER   
            The number of rows in the matrix A.   

    N       (input) INTEGER   
            The number of columns in the matrix A.   

    NB      (input) INTEGER   
            The number of leading rows and columns of A to be reduced.   

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)   
            On entry, the m by n general matrix to be reduced.   
            On exit, the first NB rows and columns of the matrix are   
            overwritten; the rest of the array is unchanged.   
            If m >= n, elements on and below the diagonal in the first NB   
              columns, with the array TAUQ, represent the orthogonal   
              matrix Q as a product of elementary reflectors; and   
              elements above the diagonal in the first NB rows, with the   
              array TAUP, represent the orthogonal matrix P as a product   
              of elementary reflectors.   
            If m < n, elements below the diagonal in the first NB   
              columns, with the array TAUQ, represent the orthogonal   
              matrix Q as a product of elementary reflectors, and   
              elements on and above the diagonal in the first NB rows,   
              with the array TAUP, represent the orthogonal matrix P as   
              a product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    D       (output) COMPLEX_16 array, dimension (NB)   
            The diagonal elements of the first NB rows and columns of   
            the reduced matrix.  D(i) = A(i,i).   

    E       (output) COMPLEX_16 array, dimension (NB)   
            The off-diagonal elements of the first NB rows and columns of   
            the reduced matrix.   

    TAUQ    (output) COMPLEX_16 array dimension (NB)   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Q. See Further Details.   

    TAUP    (output) COMPLEX_16 array, dimension (NB)   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix P. See Further Details.   

    X       (output) COMPLEX_16 array, dimension (LDX,NB)   
            The m-by-nb matrix X required to update the unreduced part   
            of A.   

    LDX     (input) INTEGER   
            The leading dimension of the array X. LDX >= M.   

    Y       (output) COMPLEX_16 array, dimension (LDY,NB)   
            The n-by-nb matrix Y required to update the unreduced part   
            of A.   

    LDY     (input) INTEGER   
            The leading dimension of the array Y. LDY >= N.   

    Further Details   
    ===============   

    The matrices Q and P are represented as products of elementary   
    reflectors:   

       Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are complex scalars, and v and u are complex vectors.   

    If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in   
    A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in   
    A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).   

    If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in   
    A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in   
    A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).   

    The elements of the vectors v and u together form the m-by-nb matrix   
    V and the nb-by-n matrix U' which are needed, with X and Y, to apply   
    the transformation to the unreduced part of the matrix, using a block   
    update of the form:  A := A - V*Y' - X*U'.   

    The contents of A on exit are illustrated by the following examples   
    with nb = 2:   

    m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):   

      (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )   
      (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )   
      (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )   
      (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )   
      (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )   
      (  v1  v2  a   a   a  )   

    where a denotes an element of the original matrix which is unchanged,   
    vi denotes an element of the vector defining H(i), and ui an element   
    of the vector defining G(i).   

    =====================================================================    */


    /* Table of constant values */
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    cuDoubleComplex c_one = MAGMA_Z_ONE;
    static int c__1 = 1;
    cuDoubleComplex c_zero = MAGMA_Z_ZERO;
    
    /* System generated locals */
    int a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset, i__2, 
	    i__3;
    /* Local variables */
    static int i__;
    cuDoubleComplex alpha;

    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d;
    --e;
    --tauq;
    --taup;

    x_dim1 = ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    dx-= 1 + lddx;

    y_dim1 = ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    dy-= 1 + lddy;

    /* Function Body */
    if (m <= 0 || n <= 0) {
	return 0;
    }

    cuDoubleComplex *f = (cuDoubleComplex *)malloc(max(n,m)*sizeof(cuDoubleComplex ));
    static cudaStream_t stream[2];
    cudaStreamCreate(&stream[0]);
    cudaStreamCreate(&stream[1]);

    if (m >= n) {

        /* Reduce to upper bidiagonal form */

	for (i__ = 1; i__ <= nb; ++i__) {

	    /*  Update A(i:m,i) */
	    i__2 = m - i__ + 1;
	    i__3 = i__ - 1;
#if defined(PRECISION_z) || defined(PRECISION_c)
            lapackf77_zlacgv( &i__3, &y[i__+y_dim1], &ldy );
#endif
	    blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one, &a[i__ + a_dim1], &lda,
		   &y[i__+y_dim1], &ldy, &c_one, &a[i__ + i__ * a_dim1], &c__1);
#if defined(PRECISION_z) || defined(PRECISION_c)
            lapackf77_zlacgv( &i__3, &y[i__+y_dim1], &ldy );
#endif
	    blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one, &x[i__ + x_dim1], &ldx, 
		   &a[i__*a_dim1+1], &c__1, &c_one, &a[i__+i__*a_dim1], &c__1);
	    
	    /* Generate reflection Q(i) to annihilate A(i+1:m,i) */

            alpha = a[i__ + i__ * a_dim1];
	    i__2 = m - i__ + 1;
	    i__3 = i__ + 1;
	    lapackf77_zlarfg(&i__2, &alpha, 
		    &a[min(i__3,m) + i__ * a_dim1], &c__1, &tauq[i__]);
	    d[i__] = MAGMA_Z_GET_X( alpha );
	    if (i__ < n) {
		a[i__ + i__ * a_dim1] = c_one;

		/* Compute Y(i+1:n,i) */
		i__2 = m - i__ + 1;
		i__3 = n - i__;

		// 1. Send the block reflector  A(i+1:m,i) to the GPU ------
		cublasSetVector(i__2, sizeof(cuDoubleComplex),
				a + i__   + i__   * a_dim1, 1,
				da+(i__-1)+(i__-1)* (ldda), 1);
		// 2. Multiply ---------------------------------------------
		cublasZgemv(MagmaConjTrans, i__2, i__3, c_one, 
			    da + (i__-1) + ((i__-1) + 1) * (ldda), ldda, 
			    da + (i__-1) + (i__-1) * (ldda), c__1, c_zero, 
			    dy + i__ + 1 + i__ * y_dim1, c__1);
		
		// 3. Put the result back ----------------------------------
		cudaMemcpy2DAsync(y+i__+1+i__*y_dim1, y_dim1*sizeof(cuDoubleComplex),
				  dy+i__+1+i__*y_dim1, y_dim1*sizeof(cuDoubleComplex),
				  sizeof(cuDoubleComplex)*i__3, 1,
				  cudaMemcpyDeviceToHost,stream[1]);
		i__2 = m - i__ + 1;
		i__3 = i__ - 1;
		blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_one, &a[i__ + a_dim1], 
			&lda, &a[i__ + i__ * a_dim1], &c__1, &c_zero, 
		       &y[i__ * y_dim1 + 1], &c__1);

		i__2 = n - i__;
                i__3 = i__ - 1;
                blasf77_zgemv("N", &i__2, &i__3, &c_neg_one, &y[i__ + 1 +y_dim1], &ldy,
		       &y[i__ * y_dim1 + 1], &c__1,
		       &c_zero, f, &c__1);
                i__2 = m - i__ + 1;
                i__3 = i__ - 1;
                blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_one, &x[i__ + x_dim1],
		       &ldx, &a[i__ + i__ * a_dim1], &c__1, &c_zero,
		       &y[i__ * y_dim1 + 1], &c__1);
		
		// 4. Synch to make sure the result is back ----------------
		cudaStreamSynchronize(stream[1]);

		if (i__3!=0){
		  i__2 = n - i__;
		  blasf77_zaxpy(&i__2, &c_one, f,&c__1, &y[i__+1+i__*y_dim1],&c__1);
		}

		i__2 = i__ - 1;
		i__3 = n - i__;
		blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_neg_one, &a[(i__ + 1) *
			a_dim1 + 1], &lda, &y[i__ * y_dim1 + 1], &c__1, &c_one,
			&y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = n - i__;
		blasf77_zscal(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);

		/* Update A(i,i+1:n) */
		i__2 = n - i__;
#if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_zlacgv( &i__2, &a[i__+(i__+1)*a_dim1], &lda );
                lapackf77_zlacgv( &i__,  &a[i__+a_dim1], &lda );
#endif
		blasf77_zgemv("No transpose", &i__2, &i__, &c_neg_one, &y[i__ + 1 +
			y_dim1], &ldy, &a[i__ + a_dim1], &lda, &c_one, &a[i__ + (
			i__ + 1) * a_dim1], &lda);
		i__2 = i__ - 1;
		i__3 = n - i__;
#if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_zlacgv( &i__,  &a[i__+a_dim1], &lda );
                lapackf77_zlacgv( &i__2, &x[i__+x_dim1], &ldx );
#endif
		blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_neg_one, &a[(i__ + 1) *
			a_dim1 + 1], &lda, &x[i__ + x_dim1], &ldx, &c_one, &a[
			i__ + (i__ + 1) * a_dim1], &lda);
#if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_zlacgv( &i__2, &x[i__+x_dim1], &ldx );
#endif

		/* Generate reflection P(i) to annihilate A(i,i+2:n) */
		i__2 = n - i__;
		/* Computing MIN */
		i__3 = i__ + 2;
		alpha = a[i__ + (i__ + 1) * a_dim1];
		lapackf77_zlarfg(&i__2, &alpha, &a[i__ + min(
			i__3,n) * a_dim1], &lda, &taup[i__]);
		e[i__] = MAGMA_Z_GET_X ( alpha );
		a[i__ + (i__ + 1) * a_dim1] = c_one;

		/* Compute X(i+1:m,i) */
		i__2 = m - i__;
		i__3 = n - i__;
                // 1. Send the block reflector  A(i+1:m,i) to the GPU ------
                cublasSetVector(i__3, sizeof(cuDoubleComplex),
                                a + i__   + (i__   +1)* a_dim1, lda,
                                da+(i__-1)+((i__-1)+1)*(ldda), ldda);
                // 2. Multiply ---------------------------------------------
		//cublasZcopy(i__3, da+(i__-1)+((i__-1)+1)*(ldda), ldda,
		//	    dy + 1 + lddy, 1);
                cublasZgemv('N', i__2, i__3, c_one,
                            da + (i__-1)+1+ ((i__-1)+1) * (ldda), ldda,
                            da + (i__-1) +  ((i__-1)+1) * (ldda), ldda,
			    //dy + 1 + lddy, 1,
			    c_zero, dx + i__ + 1 + i__ * x_dim1, c__1);

		// 3. Put the result back ----------------------------------
		cudaMemcpy2DAsync(x+i__+1+i__*x_dim1, x_dim1*sizeof(cuDoubleComplex),
				  dx+i__+1+i__*x_dim1, x_dim1*sizeof(cuDoubleComplex),
				  sizeof(cuDoubleComplex)*i__2, 1,
                                  cudaMemcpyDeviceToHost,stream[1]);

		i__2 = n - i__;
		blasf77_zgemv(MagmaConjTransStr, &i__2, &i__, &c_one, &y[i__ + 1 + y_dim1],
			&ldy, &a[i__ + (i__ + 1) * a_dim1], &lda, &c_zero, &x[
			i__ * x_dim1 + 1], &c__1);

		i__2 = m - i__;
                blasf77_zgemv("N", &i__2, &i__, &c_neg_one, &a[i__ + 1 + a_dim1], &lda,
		       &x[i__ * x_dim1 + 1], &c__1, &c_zero, f, &c__1);
                i__2 = i__ - 1;
                i__3 = n - i__;
		blasf77_zgemv("N", &i__2, &i__3, &c_one, &a[(i__ + 1) * a_dim1 + 1],
		       &lda, &a[i__ + (i__ + 1) * a_dim1], &lda,
		       &c_zero, &x[i__ * x_dim1 + 1], &c__1);

		// 4. Synch to make sure the result is back ----------------
                cudaStreamSynchronize(stream[1]);
		if (i__!=0){
                  i__2 = m - i__;
                  blasf77_zaxpy(&i__2, &c_one, f,&c__1, &x[i__+1+i__*x_dim1],&c__1);
		}


		i__2 = m - i__;
		i__3 = i__ - 1;
		blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one, &x[i__ + 1 + 
			x_dim1], &ldx, &x[i__ * x_dim1 + 1], &c__1, &c_one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = m - i__;
		blasf77_zscal(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);

#if defined(PRECISION_z) || defined(PRECISION_c)
		i__2 = n - i__;
                lapackf77_zlacgv( &i__2,  &a[i__+(i__+1)*a_dim1], &lda );
                // 4. Send the block reflector  A(i+1:m,i) to the GPU after ZLACGV()
                cublasSetVector(i__2, sizeof(cuDoubleComplex),
                                a + i__   + (i__   +1)* a_dim1, lda,
                                da+(i__-1)+((i__-1)+1)*(ldda), ldda);
#endif
	    }
	}
    } else {

      /* Reduce to lower bidiagonal form */

      for (i__ = 1; i__ <= nb; ++i__) {

	/* Update A(i,i:n) */
	i__2 = n - i__ + 1;
	i__3 = i__ - 1;
#if defined(PRECISION_z) || defined(PRECISION_c)
	lapackf77_zlacgv(&i__2, &a[i__ + i__ * a_dim1], &lda);
	lapackf77_zlacgv(&i__3, &a[i__ + a_dim1], &lda);
#endif
	blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one, &y[i__ + y_dim1], &ldy,
	       &a[i__ + a_dim1], &lda, &c_one, &a[i__ + i__ * a_dim1], &lda);
        i__2 = i__ - 1;
#if defined(PRECISION_z) || defined(PRECISION_c)
	lapackf77_zlacgv(&i__3, &a[i__ + a_dim1], &lda);
        lapackf77_zlacgv(&i__3, &x[i__ + x_dim1], &ldx);
#endif
	i__3 = n - i__ + 1;
	blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_neg_one, &a[i__ * a_dim1 + 1],
	       &lda, &x[i__ + x_dim1], &ldx, &c_one, &a[i__ + i__ * a_dim1], &lda);
#if defined(PRECISION_z) || defined(PRECISION_c)
        lapackf77_zlacgv(&i__2, &x[i__ + x_dim1], &ldx);
#endif

	/* Generate reflection P(i) to annihilate A(i,i+1:n) */
	i__2 = n - i__ + 1;
	/* Computing MIN */
	i__3 = i__ + 1;
        alpha = a[i__ + i__ * a_dim1];
	lapackf77_zlarfg(&i__2, &alpha, 
		&a[i__ + min(i__3,n) * a_dim1], &lda, &taup[i__]);
	d[i__] = MAGMA_Z_GET_X( alpha );
	if (i__ < m) {
	  a[i__ + i__ * a_dim1] = c_one;
	  
	  /* Compute X(i+1:m,i) */
	  i__2 = m - i__;
	  i__3 = n - i__ + 1;
	  blasf77_zgemv("No transpose", &i__2, &i__3, &c_one,
		 &a[i__ + 1 + i__ * a_dim1], &lda, &a[i__ + i__ * a_dim1], &lda,
		 &c_zero, &x[i__ + 1 + i__ * x_dim1], &c__1);
	  i__2 = n - i__ + 1;
	  i__3 = i__ - 1;
	  blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_one, &y[i__ + y_dim1],
		 &ldy, &a[i__ + i__ * a_dim1], &lda, &c_zero,
		 &x[i__ *  x_dim1 + 1], &c__1);
	  i__2 = m - i__;
	  i__3 = i__ - 1;
	  blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one,
		 &a[i__ + 1 + a_dim1], &lda, &x[i__ * x_dim1 + 1], &c__1, &c_one,
		 &x[i__ + 1 + i__ * x_dim1], &c__1);
	  i__2 = i__ - 1;
	  i__3 = n - i__ + 1;
	  blasf77_zgemv("No transpose", &i__2, &i__3, &c_one,
		 &a[i__ * a_dim1 + 1], &lda, &a[i__ + i__ * a_dim1], &lda, &c_zero,
		 &x[i__ * x_dim1 + 1], &c__1);
	  i__2 = m - i__;
	  i__3 = i__ - 1;
	  blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one, 
		 &x[i__ + 1 + x_dim1], &ldx, &x[i__ * x_dim1 + 1], &c__1, &c_one,
		 &x[i__ + 1 + i__ * x_dim1], &c__1);
	  i__2 = m - i__;
	  blasf77_zscal(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
	  i__2 = n - i__ + 1;
#if defined(PRECISION_z) || defined(PRECISION_c)
	  lapackf77_zlacgv(&i__2, &a[i__ + i__ * a_dim1], &lda);
#endif
	  
	  /* Update A(i+1:m,i) */
	  i__2 = m - i__;
	  i__3 = i__ - 1;
#if defined(PRECISION_z) || defined(PRECISION_c)
	  lapackf77_zlacgv(&i__3, &y[i__ + y_dim1], &ldy);
#endif
	  blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one, 
		 &a[i__ + 1 + a_dim1], &lda, &y[i__ + y_dim1], &ldy, &c_one, 
		 &a[i__ + 1 + i__ * a_dim1], &c__1);
	  i__2 = m - i__;
#if defined(PRECISION_z) || defined(PRECISION_c)
	  lapackf77_zlacgv(&i__3, &y[i__ + y_dim1], &ldy);
#endif
	  blasf77_zgemv("No transpose", &i__2, &i__, &c_neg_one, 
		 &x[i__ + 1 + x_dim1], &ldx, &a[i__ * a_dim1 + 1], &c__1, &c_one,
		 &a[i__ + 1 + i__ * a_dim1], &c__1);
	  
	  /* Generate reflection Q(i) to annihilate A(i+2:m,i) */
	  i__2 = m - i__;
	  /* Computing MIN */
	  i__3 = i__ + 2;
          alpha = a[i__ + 1 + i__ * a_dim1];
	  lapackf77_zlarfg(&i__2, &alpha,
		  &a[min(i__3,m) + i__ * a_dim1], &c__1, &tauq[i__]);
	  e[i__] = MAGMA_Z_GET_X( alpha );
	  a[i__ + 1 + i__ * a_dim1] = c_one;
	  
	  /* Compute Y(i+1:n,i) */
	  i__2 = m - i__;
	  i__3 = n - i__;
	  blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_one, 
		 &a[i__ + 1 + (i__ + 1) * a_dim1], &lda, 
		 &a[i__ + 1 + i__ * a_dim1], &c__1, 
		 &c_zero, &y[i__ + 1 + i__ * y_dim1], &c__1);
	  i__2 = m - i__;
	  i__3 = i__ - 1;
	  blasf77_zgemv(MagmaConjTransStr, &i__2, &i__3, &c_one, &a[i__ + 1 + a_dim1], 
		 &lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_zero, 
		 &y[ i__ * y_dim1 + 1], &c__1);
	  i__2 = n - i__;
	  i__3 = i__ - 1;
	  blasf77_zgemv("No transpose", &i__2, &i__3, &c_neg_one,
		 &y[i__ + 1 + y_dim1], &ldy, &y[i__ * y_dim1 + 1], &c__1,
		 &c_one, &y[i__ + 1 + i__ * y_dim1], &c__1);
	  i__2 = m - i__;
	  blasf77_zgemv(MagmaConjTransStr, &i__2, &i__, &c_one, &x[i__ + 1 + x_dim1],
		 &ldx, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_zero,
		 &y[i__ * y_dim1 + 1], &c__1);
	  i__2 = n - i__;
	  blasf77_zgemv(MagmaConjTransStr, &i__, &i__2, &c_neg_one,
		 &a[(i__ + 1) * a_dim1 + 1], &lda, &y[i__ * y_dim1 + 1],
		 &c__1, &c_one, &y[i__ + 1 + i__ * y_dim1], &c__1);
	  i__2 = n - i__;
	  blasf77_zscal(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
#if defined(PRECISION_z) || defined(PRECISION_c)
	} else {
          i__2 = n - i__ + 1;
          lapackf77_zlacgv(&i__2, &a[i__ + i__ * a_dim1], &lda);
#endif
	}
      }
    }
    
    free(f);
    
    return 0;

    /* End of MAGMA_ZLABRD */

} /* zlabrd_ */

