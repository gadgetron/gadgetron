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
magma_zungqr(magma_int_t m, magma_int_t n, magma_int_t k,
	     cuDoubleComplex *a, magma_int_t lda,
	     cuDoubleComplex *tau, cuDoubleComplex *dT,
	     magma_int_t nb, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    ZUNGQR generates an M-by-N COMPLEX_16 matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by ZGEQRF.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix Q. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix Q. M >= N >= 0.

    K       (input) INTEGER
            The number of elementary reflectors whose product defines the
            matrix Q. N >= K >= 0.

    NB      (input) INTEGER
            The block size used in the generation of the elementary
            reflectors H(i) in DA. Thus this gives the dimensions of
            the matrices T, stored in DT.

    A       (input/output) COMPLEX_16 array A, dimension (LDDA,N). 
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by ZGEQRF_GPU in the 
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    LDA     (input) INTEGER
            The first dimension of the array A. LDA >= max(1,M).

    TAU     (input) COMPLEX_16 array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZGEQRF_GPU.

    DT      (input) COMPLEX_16 array on the GPU device.
            DT contains the T matrices used in blocking the elementary
            reflectors H(i), e.g., this can be the 6th argument of 
            magma_zgeqrf_gpu.

    NB      (input) INTEGER
            This is the block size used in ZGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument has an illegal value
    =====================================================================    */

    #define  a_ref(i,j)     ( a + (j)*lda  + (i))
    #define da_ref(i,j)     (da + (j)*ldda + (i))
    #define t_ref(a_1)      (dT+(a_1)*nb)

    magma_int_t  i__1, i__2, i__3;
    magma_int_t lwork, ldda;
    static magma_int_t i, ib, ki, kk, iinfo;
    magma_int_t lddwork = min(m, n);
    cuDoubleComplex *da, *work, *dwork;
    static cudaStream_t stream;

    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if ((n < 0) || (n > m)) {
	*info = -2;
    } else if ((k < 0) || (k > n)) {
	*info = -3;
    } else if (lda < max(1,m)) {
	*info = -5;
    }
    if (*info != 0)
      return MAGMA_ERR_ILLEGAL_VALUE;

    if (n <= 0)
      return MAGMA_SUCCESS;

    /* Allocate GPU work space */
    ldda = ((m+31)/32)*32;
    lddwork = ((lddwork+31)/32)*32;
    if (CUBLAS_STATUS_SUCCESS != 
	cublasAlloc((n)*ldda + nb*lddwork, sizeof(cuDoubleComplex), (void**)&da)) 
      {
	*info = -11;
	return MAGMA_ERR_CUBLASALLOC;
      }
    dwork = da + (n)*ldda;

    /* Allocate CPU work space */
    lwork = n * nb;
    work = (cuDoubleComplex *)malloc(lwork*sizeof(cuDoubleComplex));
    if( work == NULL ) {
        cublasFree(da);
        magma_xerbla(__func__, info);
	return MAGMA_ERR_ALLOCATION;
    }

    cudaStreamCreate(&stream);

    if ( (nb > 1) && (nb < k) )
      {
	/*  Use blocked code after the last block.
	    The first kk columns are handled by the block method. */
	ki = (k - nb - 1) / nb * nb;
	kk = min(k, ki + nb);

	/* Set A(1:kk,kk+1:n) to zero. */
        magmablas_zlaset(kk, n-kk, da_ref(0,kk), ldda);
      }
    else
      kk = 0;

    /* Use unblocked code for the last or only block. */
    if (kk < n)
      {
	i__1 = m - kk;
	i__2 = n - kk;
	i__3 = k - kk;
	lapackf77_zungqr(&i__1, &i__2, &i__3, 
			 a_ref(kk, kk), &lda,
			 &tau[kk], work, &lwork, &iinfo);

	cublasSetMatrix(i__1, i__2, sizeof(cuDoubleComplex),
 			 a_ref(kk, kk), lda, 
			da_ref(kk, kk), ldda);
      }

    if (kk > 0)
      {
	/* Use blocked code */
	for (i = ki; i >= 0; i-=nb)
	  {
	    ib = min(nb, k - i);

	    /* Send the current panel to the GPU */
	    i__2 = m - i;
	    zpanel_to_q(MagmaUpper, ib, a_ref(i,i), lda, work);
	    cublasSetMatrix( i__2, ib, sizeof(cuDoubleComplex),
			      a_ref(i, i), lda,
			     da_ref(i, i), ldda);
			     
	    if (i + ib < n)
	      {
		/* Apply H to A(i:m,i+ib:n) from the left */
		i__3 = n - i - ib;
		magma_zlarfb_gpu( MagmaLeft, MagmaNoTrans, MagmaForward, MagmaColumnwise,
				  i__2, i__3, ib,
				  da_ref(i, i   ), ldda, t_ref(i),      nb,
				  da_ref(i, i+ib), ldda,    dwork, lddwork);
	      }

	    /* Apply H to rows i:m of current block on the CPU */
	    lapackf77_zungqr(&i__2, &ib, &ib, 
			     a_ref(i, i), &lda, 
			     &tau[i], work, &lwork, &iinfo);
	    cudaMemcpy2DAsync(da_ref(i,i), ldda * sizeof(cuDoubleComplex),
			       a_ref(i,i), lda * sizeof(cuDoubleComplex),
			      sizeof(cuDoubleComplex)*i__2, ib,
                              cudaMemcpyHostToDevice, stream);

	    /* Set rows 1:i-1 of current block to zero */
            i__2 = i + ib;
	    magmablas_zlaset(i, i__2 - i, da_ref(0,i), ldda);
	  }
      }

    cublasGetMatrix(m, n, sizeof(cuDoubleComplex),
		    da_ref(0, 0), ldda, a_ref(0, 0), lda);

    
    cudaStreamDestroy(stream);
    cublasFree(da);
    free(work);

    return MAGMA_SUCCESS;
} /* magma_zungqr */

#undef da_ref
#undef a_ref
#undef t_ref
