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
magma_ztsqrt_gpu(int *m, int *n, 
		 cuDoubleComplex *a1, cuDoubleComplex *a2, int  *lda,
		 cuDoubleComplex  *tau, cuDoubleComplex *work, 
		 int *lwork, cuDoubleComplex *dwork, int *info )
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   

    ZGEQRF computes a QR factorization of a complex M-by-N matrix A:   
    A = Q * R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX_16 array on the GPU, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the elements on and above the diagonal of the array   
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is   
            upper triangular if m >= n); the elements below the diagonal,   
            with the array TAU, represent the orthogonal matrix Q as a   
            product of min(m,n) elementary reflectors (see Further   
            Details).

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    TAU     (output) COMPLEX_16 array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    WORK    (workspace/output) COMPLEX_16 array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

            Higher performance is achieved if WORK is in pinned memory, e.g.
            allocated using cudaMallocHost.

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= (M+N+NB)*NB,   
            where NB can be obtained through magma_get_zgeqrf_nb(M).

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued.   

    DWORK   (workspace/output)  COMPLEX_16 array on the GPU, dimension 2*N*NB,
            where NB can be obtained through magma_get_zgeqrf_nb(M).
            It starts with NB*NB blocks that store the triangular T 
            matrices, followed by the NB*NB blocks of the diagonal 
            inverses for the R matrix.

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),   
    and tau in TAU(i).   
    =====================================================================    */

   #define a1_ref(a_1,a_2) ( a1+(a_2)*(*lda) + (a_1))
   #define a2_ref(a_1,a_2) ( a2+(a_2)*(*lda) + (a_1))
   #define t_ref(a_1)     (dwork+(a_1))
   #define d_ref(a_1)     (dwork+(lddwork+(a_1))*nb)
   #define dd_ref(a_1)    (dwork+(2*lddwork+(a_1))*nb)
   #define work_a1        ( work )
   #define work_a2        ( work + nb ) 
   #define hwork          ( work + (nb)*(*m))

   int i, k, ldwork, lddwork, old_i, old_ib, rows, cols;
   int nbmin, ib, ldda;

   /* Function Body */
   *info = 0;
   int nb = magma_get_zgeqrf_nb(*m);

   int lwkopt = (*n+*m) * nb;
   work[0] = (cuDoubleComplex) lwkopt;
   long int lquery = *lwork == -1;
   if (*m < 0) {
     *info = -1;
   } else if (*n < 0) {
     *info = -2;
   } else if (*lda < max(1,*m)) {
     *info = -4;
   } else if (*lwork < max(1,*n) && ! lquery) {
     *info = -7;
   }
   if (*info != 0)
     return 0;
   else if (lquery)
     return 0;

   k = min(*m,*n);
   if (k == 0) {
     work[0] = 1.f;
     return 0;
   }

   int lhwork = *lwork - (*m)*nb;

   static cudaStream_t stream[2];
   cudaStreamCreate(&stream[0]);
   cudaStreamCreate(&stream[1]);

   ldda = *m;
   nbmin = 2;
   ldwork = *m;
   lddwork= k;

   // This is only blocked code for now
   for (i = 0; i < k; i += nb) {
     ib = min(k-i, nb);
     rows = *m -i;
     rows = *m;
     // Send the next panel (diagonal block of A1 & block column of A2) 
     // to the CPU (in work_a1 and work_a2)
     cudaMemcpy2DAsync(  work_a2, ldwork*sizeof(cuDoubleComplex),
			 a2_ref(0,i), (*lda)*sizeof(cuDoubleComplex),
			 sizeof(cuDoubleComplex)*rows, ib,
			 cudaMemcpyDeviceToHost,stream[1]);
     cudaMemcpy2DAsync(  work_a1, ldwork*sizeof(cuDoubleComplex),
			 // a1_ref(i,i), (*lda)*sizeof(cuDoubleComplex),
			 // the diagonal of a1 is in d_ref generated and
			 // passed from magma_zgeqrf_gpu
			 d_ref(i), ib*sizeof(cuDoubleComplex),
			 sizeof(cuDoubleComplex)*ib, ib,
			 cudaMemcpyDeviceToHost,stream[1]);
     
	if (i>0) {
	  /* Apply H' to A(i:m,i+2*ib:n) from the left */
	  // update T2
	  cols = *n-old_i-2*old_ib;
	  magma_zssrfb(*m, cols, &old_ib,
		       a2_ref(    0, old_i), lda, t_ref(old_i), &lddwork,
		       a1_ref(old_i, old_i+2*old_ib), lda, 
                       a2_ref(    0, old_i+2*old_ib), lda,
                       dd_ref(0), &lddwork);
	}

	cudaStreamSynchronize(stream[1]);

	// TTT - here goes the CPU PLASMA code
	//       Matrix T has to be put in hwork with lda = ib and 0s
	//       in the parts that are not used - copied on GPU in t_ref(i) 
	
	// Now diag of A1 is updated, send it back asynchronously to the GPU.
        // We have to play interchaning these copies to see which is faster
	cudaMemcpy2DAsync(d_ref(i), ib * sizeof(cuDoubleComplex),
			  work_a1 , ib * sizeof(cuDoubleComplex),
			  sizeof(cuDoubleComplex)*ib, ib,
			  cudaMemcpyHostToDevice,stream[0]);
	// Send the panel from A2 back to the GPU
	cublasSetMatrix(*m, ib, sizeof(cuDoubleComplex),
			work_a2, ldwork, a2_ref(0,i), *lda);

	if (i + ib < *n) {
	  // Send the triangular factor T from hwork to the GPU in t_ref(i)
	  cublasSetMatrix(ib, ib, sizeof(cuDoubleComplex), hwork, ib, t_ref(i), lddwork);

	  if (i+nb < k){
	    /* Apply H' to A(i:m,i+ib:i+2*ib) from the left */
	    // if we can do one more step, first update T1
	    magma_zssrfb(*m, ib, &ib, 
			 a2_ref(    0, i), lda, t_ref(i), &lddwork,
			 a1_ref(    i, i+ib), lda,
			 a2_ref(    0, i+ib), lda,
			 dd_ref(0), &lddwork);
	  }
	  else {
	    cols = *n-i-ib;
	    // otherwise, update until the end and fix the panel
	    magma_zssrfb(*m, cols, &ib, 
			 a2_ref(    0, i), lda, t_ref(i), &lddwork,
			 a1_ref(    i, i+ib), lda,
                         a2_ref(    0, i+ib), lda,
                         dd_ref(0), &lddwork);
	  }
	  old_i = i;
	  old_ib = ib;
	}
   }  
   
   return 0; 
   
   /* End of MAGMA_ZTSQRT_GPU */

} /* magma_ztsqrt_gpu */

#undef a1_ref
#undef a2_ref
#undef t_ref
#undef d_ref
#undef dd_ref
#undef hwork 
#undef work_a1
#undef work_a2
