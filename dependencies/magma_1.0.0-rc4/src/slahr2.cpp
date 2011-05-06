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
// #define cublasSgemv magmablas_sgemv
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t 
magma_slahr2(magma_int_t n, magma_int_t k, magma_int_t nb,
	     float *da, float *dv, 
	     float *a, magma_int_t lda,
	     float *tau, float *t, magma_int_t ldt, 
	     float *y, magma_int_t ldy)
{
/*  -- MAGMA auxiliary routine (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   

    SLAHR2 reduces the first NB columns of a real general n-BY-(n-k+1)   
    matrix A so that elements below the k-th subdiagonal are zero. The   
    reduction is performed by an orthogonal similarity transformation   
    Q' * A * Q. The routine returns the matrices V and T which determine   
    Q as a block reflector I - V*T*V', and also the matrix Y = A * V.   

    This is an auxiliary routine called by SGEHRD.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.   

    K       (input) INTEGER   
            The offset for the reduction. Elements below the k-th   
            subdiagonal in the first NB columns are reduced to zero.   
            K < N.   

    NB      (input) INTEGER   
            The number of columns to be reduced.

    DA      (input/output) REAL array on the GPU, dimension (LDA,N-K+1)   
            On entry, the n-by-(n-k+1) general matrix A.   
            On exit, the elements on and above the k-th subdiagonal in   
            the first NB columns are overwritten with the corresponding   
            elements of the reduced matrix; the elements below the k-th   
            subdiagonal, with the array TAU, represent the matrix Q as a   
            product of elementary reflectors. The other columns of A are   
            unchanged. See Further Details.   

    DV      (output) REAL array on the GPU, dimension (N, NB)
            On exit this contains the Householder vectors of the transformation.

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    TAU     (output) REAL array, dimension (NB)   
            The scalar factors of the elementary reflectors. See Further   
            Details.   

    T       (output) REAL array, dimension (LDT,NB)   
            The upper triangular matrix T.   

    LDT     (input) INTEGER   
            The leading dimension of the array T.  LDT >= NB.   

    Y       (output) REAL array, dimension (LDY,NB)   
            The n-by-nb matrix Y.   

    LDY     (input) INTEGER   
            The leading dimension of the array Y. LDY >= N.   

    Further Details   
    ===============   
    The matrix Q is represented as a product of nb elementary reflectors   

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in   
    A(i+k+1:n,i), and tau in TAU(i).   

    The elements of the vectors v together form the (n-k+1)-by-nb matrix   
    V which is needed, with T and Y, to apply the transformation to the   
    unreduced part of the matrix, using an update of the form:   
    A := (I - V*T*V') * (A - Y*T*V').   

    The contents of A on exit are illustrated by the following example   
    with n = 7, k = 3 and nb = 2:   

       ( a   a   a   a   a )   
       ( a   a   a   a   a )   
       ( a   a   a   a   a )   
       ( h   h   a   a   a )   
       ( v1  h   a   a   a )   
       ( v1  v2  a   a   a )   
       ( v1  v2  a   a   a )   

    where a denotes an element of the original matrix A, h denotes a   
    modified element of the upper Hessenberg matrix H, and vi denotes an   
    element of the vector defining H(i).

    This implementation follows the hybrid algorithm and notations described in

    S. Tomov and J. Dongarra, "Accelerating the reduction to upper Hessenberg
    form through hybrid GPU-based computing," University of Tennessee Computer
    Science Technical Report, UT-CS-09-642 (also LAPACK Working Note 219),
    May 24, 2009.
    =====================================================================    */


    float c_zero    = MAGMA_S_ZERO;
    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;

    magma_int_t ldda = n;
    magma_int_t c__1 = 1;
    
    magma_int_t a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__2, i__3;
    float d__1;

    magma_int_t i__;
    float ei;

    --tau;
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    y_dim1 = ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    if (n <= 1)
      return 0;
    
    for (i__ = 1; i__ <= nb; ++i__) {
	if (i__ > 1) {

	  /* Update A(K+1:N,I); Update I-th column of A - Y * V' */
	  i__2 = n - k + 1;
	  i__3 = i__ - 1;
          #if defined(PRECISION_z) || defined(PRECISION_c)
	     lapackf77_slacgv(&i__3, &a[k+i__-1+a_dim1], &lda);
          #endif
	  blasf77_scopy(&i__3, &a[k+i__-1+a_dim1], &lda, &t[nb*t_dim1+1], &c__1);
	  blasf77_strmv("u","n","n",&i__3,&t[t_offset], &ldt, &t[nb*t_dim1+1], &c__1);

	  blasf77_sgemv("NO TRANSPOSE", &i__2, &i__3, &c_neg_one, &y[k + y_dim1],
			&ldy, &t[nb*t_dim1+1], &c__1, &c_one, &a[k+i__*a_dim1],&c__1);

          #if defined(PRECISION_z) || defined(PRECISION_c)
             lapackf77_slacgv(&i__3, &a[k+i__-1+a_dim1], &lda);
          #endif

	  /* Apply I - V * T' * V' to this column (call it b) from the   
             left, using the last column of T as workspace   

             Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)   
                      ( V2 )             ( b2 )   
             where V1 is unit lower triangular   
             w := V1' * b1                                                 */
	  
	  i__2 = i__ - 1;
	  blasf77_scopy(&i__2, &a[k+1+i__*a_dim1], &c__1, &t[nb*t_dim1+1], &c__1);
	  blasf77_strmv("Lower", MagmaTransStr, "UNIT", &i__2, 
			&a[k + 1 + a_dim1], &lda, &t[nb * t_dim1 + 1], &c__1);

	  /* w := w + V2'*b2 */
	  i__2 = n - k - i__ + 1;
	  i__3 = i__ - 1;
	  blasf77_sgemv(MagmaTransStr, &i__2, &i__3, &c_one, 
			&a[k + i__ + a_dim1], &lda, &a[k+i__+i__*a_dim1], &c__1, 
			&c_one, &t[nb*t_dim1+1], &c__1);

	  /* w := T'*w */
	  i__2 = i__ - 1;
	  blasf77_strmv("U", MagmaTransStr, "N", &i__2, &t[t_offset], &ldt, 
			&t[nb*t_dim1+1], &c__1);
	  
	  /* b2 := b2 - V2*w */
	  i__2 = n - k - i__ + 1;
	  i__3 = i__ - 1;
	  blasf77_sgemv("N", &i__2, &i__3, &c_neg_one, &a[k + i__ + a_dim1], &lda, 
		 &t[nb*t_dim1+1], &c__1, &c_one, &a[k+i__+i__*a_dim1], &c__1);

	  /* b1 := b1 - V1*w */
	  i__2 = i__ - 1;
	  blasf77_strmv("L","N","U",&i__2,&a[k+1+a_dim1],&lda,&t[nb*t_dim1+1],&c__1);
	  blasf77_saxpy(&i__2, &c_neg_one, &t[nb * t_dim1 + 1], &c__1, 
		 &a[k + 1 + i__ * a_dim1], &c__1);
	  
	  a[k + i__ - 1 + (i__ - 1) * a_dim1] = ei;
	}
	
	/* Generate the elementary reflector H(I) to annihilate A(K+I+1:N,I) */
	i__2 = n - k - i__ + 1;
	i__3 = k + i__ + 1;
	lapackf77_slarfg(&i__2, &a[k + i__ + i__ * a_dim1], 
                         &a[min(i__3,n) + i__ * a_dim1], &c__1, &tau[i__]);
	ei = a[k + i__ + i__ * a_dim1];
	a[k + i__ + i__ * a_dim1] = c_one;

	/* Compute  Y(K+1:N,I) */
        i__2 = n - k;
	i__3 = n - k - i__ + 1;
        cublasSetVector(i__3, sizeof(float), 
                        &a[k + i__ + i__*a_dim1], 1, dv+(i__-1)*(ldda+1), 1);

	cublasSgemv('N', i__2+1, i__3, c_one, 
		    da -1 + k + i__ * ldda, ldda, 
		    dv+(i__-1)*(ldda+1), c__1, c_zero, 
		    da-1 + k + (i__-1)*ldda, c__1);     
	
	i__2 = n - k - i__ + 1;
	i__3 = i__ - 1;
	blasf77_sgemv(MagmaTransStr, &i__2, &i__3, &c_one, 
		      &a[k + i__ + a_dim1], &lda, &a[k+i__+i__*a_dim1], &c__1, 
		      &c_zero, &t[i__*t_dim1+1], &c__1);

	/* Compute T(1:I,I) */
	i__2 = i__ - 1;
	MAGMA_S_OP_NEG_ASGN( d__1, tau[i__] );
	blasf77_sscal(&i__2, &d__1, &t[i__ * t_dim1 + 1], &c__1);
	blasf77_strmv("U","N","N", &i__2, &t[t_offset], &ldt, &t[i__*t_dim1+1], &c__1);
	t[i__ + i__ * t_dim1] = tau[i__];

        cublasGetVector(n - k + 1, sizeof(float),
	                da-1+ k+(i__-1)*ldda, 1, y+ k + i__*y_dim1, 1);
    }
    a[k + nb + nb * a_dim1] = ei;

    return 0;
} /* magma_slahr2 */

