/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d

*/
#include "common_magma.h"

#define PRECISION_d

extern "C" magma_int_t 
magma_dlatrd(char uplo, magma_int_t n_, magma_int_t nb_, 
	     double *a,  magma_int_t lda_, 
	     double *e, double *tau, 
	     double *w,  magma_int_t ldw_,
	     double *da, magma_int_t ldda_, 
	     double *dw, magma_int_t lddw_)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   
    DLATRD reduces NB rows and columns of a real symmetric matrix A to   
    symmetric tridiagonal form by an orthogonal similarity   
    transformation Q' * A * Q, and returns the matrices V and W which are   
    needed to apply the transformation to the unreduced part of A.   

    If UPLO = 'U', DLATRD reduces the last NB rows and columns of a   
    matrix, of which the upper triangle is supplied;   
    if UPLO = 'L', DLATRD reduces the first NB rows and columns of a   
    matrix, of which the lower triangle is supplied.   

    This is an auxiliary routine called by DSYTRD.   

    Arguments   
    =========   
    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U': Upper triangular   
            = 'L': Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.   

    NB      (input) INTEGER   
            The number of rows and columns to be reduced.   

    A       (input/output) DOUBLE_PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit:   
            if UPLO = 'U', the last NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements above the diagonal   
              with the array TAU, represent the orthogonal matrix Q as a   
              product of elementary reflectors;   
            if UPLO = 'L', the first NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements below the diagonal   
              with the array TAU, represent the  orthogonal matrix Q as a   
              product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= (1,N).   

    E       (output) DOUBLE_PRECISION array, dimension (N-1)   
            If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal   
            elements of the last NB columns of the reduced matrix;   
            if UPLO = 'L', E(1:nb) contains the subdiagonal elements of   
            the first NB columns of the reduced matrix.   

    TAU     (output) DOUBLE_PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors, store     
d in   
            TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.   
            See Further Details.   

    W       (output) DOUBLE_PRECISION array, dimension (LDW,NB)   
            The n-by-nb matrix W required to update the unreduced part   
            of A.   

    LDW     (input) INTEGER   
            The leading dimension of the array W. LDW >= max(1,N).   

    Further Details   
    ===============   
    If UPLO = 'U', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(n) H(n-1) . . . H(n-nb+1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),   
    and tau in TAU(i-1).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),   
    and tau in TAU(i).   

    The elements of the vectors v together form the n-by-nb matrix V   
    which is needed, with W, to apply the transformation to the unreduced   
    part of the matrix, using a symmetric rank-2k update of the form:   
    A := A - V*W' - W*V'.   

    The contents of A on exit are illustrated by the following examples   
    with n = 5 and nb = 2:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  a   a   a   v4  v5 )              (  d                  )   
      (      a   a   v4  v5 )              (  1   d              )   
      (          a   1   v5 )              (  v1  1   a          )   
      (              d   1  )              (  v1  v2  a   a      )   
      (                  d  )              (  v1  v2  a   a   a  )   

    where d denotes a diagonal element of the reduced matrix, a denotes   
    an element of the original matrix that is unchanged, and vi denotes   
    an element of the vector defining H(i).   
    =====================================================================    */
 
  
    char uplo_[2]  = {uplo, 0};
    magma_int_t *n    = &n_;
    magma_int_t *nb   = &nb_;
    magma_int_t *lda  = &lda_;
    magma_int_t *ldw  = &ldw_;
    magma_int_t *ldda = &ldda_;
    magma_int_t *lddw = &lddw_;

    double c_neg_one = MAGMA_D_NEG_ONE;
    double c_one     = MAGMA_D_ONE;
    double c_zero    = MAGMA_D_ZERO;

    #if defined(PRECISION_z) || defined(PRECISION_c)
       double value = MAGMA_D_ZERO;
    #endif
    
    static magma_int_t c__1 = 1;

    magma_int_t a_dim1, a_offset, w_dim1, w_offset, i__2, i__3;
    static magma_int_t i__, iw;
  
    static double alpha;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --tau;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    dw-= 1 + *lddw;

    double *f = (double *)malloc((*n)*sizeof(double ));

    if (*n <= 0) {
      return 0;
    }

    static cudaStream_t stream;
    cudaStreamCreate(&stream);

    if (lapackf77_lsame(uplo_, "U")) {

      fprintf(stderr, "dlatrd: uplo has to be 'L'; 'U' not implemented \n");
      exit(-1);

      /* Reduce last NB columns of upper triangle */
      for (i__ = *n; i__ >= *n - *nb + 1; --i__) {
	iw = i__ - *n + *nb;
	if (i__ < *n) {
	  /* Update A(1:i,i) */
	  i__2 = *n - i__;
	  blasf77_dgemv("No transpose", &i__, &i__2, &c_neg_one, 
		 &a[(i__+1)*a_dim1 + 1], lda, &w[i__ + (iw + 1)*w_dim1], ldw, 
		 &c_one, &a[i__ * a_dim1 + 1], &c__1);
	  i__2 = *n - i__;
	  blasf77_dgemv("No transpose", &i__, &i__2, &c_neg_one, 
		 &w[(iw+1)*w_dim1 + 1], ldw, &a[i__ + (i__+1) * a_dim1], lda, 
		 &c_one, &a[i__ * a_dim1 + 1], &c__1);
	}
	if (i__ > 1) {
	  /* Generate elementary reflector H(i) to annihilate A(1:i-2,i) */
	  i__2 = i__ - 1;
	  lapackf77_dlarfg(&i__2, &a[i__ - 1 + i__ * a_dim1], &a[i__ * a_dim1 + 1], 
		  &c__1, &tau[i__ - 1]);

	  e[i__ - 1] = MAGMA_D_GET_X( a[i__ - 1 + i__ * a_dim1] );
	  a[i__ - 1 + i__ * a_dim1] = c_one;
  
	  /* Compute W(1:i-1,i) */
	  i__2 = i__ - 1;
	  blasf77_dsymv("Upper", &i__2, &c_one, &a[a_offset], lda, 
		 &a[i__*a_dim1 +1], &c__1, &c_zero, &w[iw* w_dim1+1], &c__1);
	  if (i__ < *n) {
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    blasf77_dgemv(MagmaTransStr, &i__2, &i__3, &c_one, 
		   &w[(iw+1)*w_dim1 + 1], ldw, &a[i__ * a_dim1 + 1], &c__1, 
		   &c_zero, &w[i__ + 1 + iw * w_dim1], &c__1);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    blasf77_dgemv("No transpose", &i__2, &i__3, &c_neg_one, 
		   &a[(i__+1)*a_dim1 + 1], lda, &w[i__ + 1 + iw * w_dim1], &
		   c__1, &c_one, &w[iw * w_dim1 + 1], &c__1);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    blasf77_dgemv(MagmaTransStr, &i__2, &i__3, &c_one, 
		   &a[(i__ + 1) * a_dim1 + 1], lda, &a[i__ * a_dim1 + 1], 
		   &c__1, &c_zero, &w[i__ + 1 + iw * w_dim1], &c__1);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    blasf77_dgemv("No transpose", &i__2, &i__3, &c_neg_one, 
		   &w[(iw + 1) *  w_dim1 + 1], ldw, &w[i__ + 1 + iw * w_dim1],
		   &c__1, &c_one, &w[iw * w_dim1 + 1], &c__1);
	  }
	  i__2 = i__ - 1;
	  blasf77_dscal(&i__2, &tau[i__ - 1], &w[iw * w_dim1 + 1], &c__1);
	  i__2 = i__ - 1;
#if defined(PRECISION_z) || defined(PRECISION_c)
	  blasf77_ddot(&value, &i__2, &w[iw*w_dim1+1], &c__1, &a[i__ * a_dim1 + 1], &c__1);
	  alpha = tau[i__ - 1] * -.5f * value;
#else
	  alpha = tau[i__ - 1] * -.5f * 
	    blasf77_ddot(&i__2, &w[iw*w_dim1+1], &c__1, &a[i__ * a_dim1 + 1], &c__1);
#endif
	  i__2 = i__ - 1;
	  blasf77_daxpy(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, 
		 &w[iw * w_dim1 + 1], &c__1);
	}
      }

    } else {

      /*  Reduce first NB columns of lower triangle */
      for (i__ = 1; i__ <= *nb; ++i__) 
	{

	  /* Update A(i:n,i) */
	  i__2 = *n - i__ + 1;
	  i__3 = i__ - 1;
          #if defined(PRECISION_z) || defined(PRECISION_c)
	      lapackf77_dlacgv(&i__3, &w[i__ + w_dim1], ldw);
          #endif
	  blasf77_dgemv("No transpose", &i__2, &i__3, &c_neg_one, 
			&a[i__ + a_dim1], lda, 
			&w[i__ + w_dim1], ldw, &c_one, 
			&a[i__ + i__ * a_dim1], &c__1);
          #if defined(PRECISION_z) || defined(PRECISION_c)
	      lapackf77_dlacgv(&i__3, &w[i__ + w_dim1], ldw);
	      lapackf77_dlacgv(&i__3, &a[i__ + a_dim1], lda);
          #endif
	  blasf77_dgemv("No transpose", &i__2, &i__3, &c_neg_one, 
			&w[i__ + w_dim1], ldw, 
			&a[i__ + a_dim1], lda, &c_one,
			&a[i__ + i__ * a_dim1], &c__1);
          #if defined(PRECISION_z) || defined(PRECISION_c)
	      lapackf77_dlacgv(&i__3, &a[i__ + a_dim1], lda);
          #endif

	  if (i__ < *n) 
	    {
	      /* Generate elementary reflector H(i) to annihilate A(i+2:n,i) */
	      i__2 = *n - i__;
	      i__3 = i__ + 2;
	      MAGMA_D_ASSIGN(alpha, a[i__ + 1 + i__ * a_dim1]);
	      lapackf77_dlarfg(&i__2, &alpha, 
			       &a[min(i__3,*n) + i__ * a_dim1], &c__1, &tau[i__]);
	      e[i__] = MAGMA_D_GET_X( alpha );
	      MAGMA_D_SET2REAL(a[i__ + 1 + i__ * a_dim1], 1.);

	      /* Compute W(i+1:n,i) */ 
	      // 1. Send the block reflector  A(i+1:n,i) to the GPU
	      cublasSetVector(i__2, sizeof(double),
			      a + i__   + 1 + i__   * a_dim1, 1,
			      da+(i__-1)+ 1 +(i__-1)* (*ldda), 1);	  
	  
	      cublasDsymv('L', i__2, c_one, da+ (i__-1)+1 + ((i__-1)+1) * (*ldda),
			  *ldda, da+ (i__-1)+1 + (i__-1)* a_dim1, c__1, c_zero,
			  dw+ i__ + 1 + i__ * w_dim1, c__1);
	  
	      // 2. Start putting the result back (asynchronously)
	      cudaMemcpy2DAsync(w +i__+1+i__*w_dim1, w_dim1*sizeof(double),
				dw+i__+1+i__*w_dim1, w_dim1*sizeof(double),
				sizeof(double)*i__2, 1,
				cudaMemcpyDeviceToHost,stream);

	      i__3 = i__ - 1;
	      blasf77_dgemv(MagmaTransStr, &i__2, &i__3, &c_one,
			    &w[i__ + 1 + w_dim1], ldw, 
			    &a[i__ + 1 + i__ * a_dim1], &c__1, &c_zero, 
			    &w[i__ * w_dim1 + 1], &c__1);

	      blasf77_dgemv("No transpose", &i__2, &i__3, &c_neg_one,
			    &a[i__ + 1 + a_dim1], lda, 
			    &w[i__ * w_dim1 + 1], &c__1,
			    &c_zero, f, &c__1);
	      blasf77_dgemv(MagmaTransStr, &i__2, &i__3, &c_one,
			    &a[i__ + 1 + a_dim1], lda, 
			    &a[i__ + 1 + i__ * a_dim1], &c__1, &c_zero,
			    &w[i__ * w_dim1 + 1], &c__1);

	      // 3. Here is where we need it
	      cudaStreamSynchronize(stream);

	      if (i__3!=0)
		blasf77_daxpy(&i__2, &c_one, f, &c__1, 
			      &w[i__ + 1 + i__ * w_dim1], &c__1);
     
	      blasf77_dgemv("No transpose", &i__2, &i__3, &c_neg_one,
			    &w[i__ + 1 + w_dim1], ldw, 
			    &w[i__ * w_dim1 + 1], &c__1, &c_one, 
			    &w[i__ + 1 + i__ * w_dim1], &c__1);
	      blasf77_dscal(&i__2, &tau[i__], &w[i__ + 1 + i__ * w_dim1], &c__1);
              #if defined(PRECISION_z) || defined(PRECISION_c)
  	          blasf77_ddot(&value, &i__2, &w[i__ +1+ i__ * w_dim1], 
				  &c__1, &a[i__ +1+ i__ * a_dim1], &c__1);
		  alpha = tau[i__]* -.5f * value;
              #else
	          alpha = tau[i__]* -.5f*blasf77_ddot(&i__2, &w[i__ +1+ i__ * w_dim1], 
						     &c__1, &a[i__ +1+ i__ * a_dim1],
						     &c__1);
              #endif
	      blasf77_daxpy(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, 
			    &w[i__ + 1 + i__ * w_dim1], &c__1);
	    }
	}
    }

    free(f);
    cudaStreamDestroy(stream);

    return 0;
} /* dlatrd_ */

