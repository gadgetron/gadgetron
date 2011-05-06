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
magma_zungqr_gpu(magma_int_t m, magma_int_t n, magma_int_t k,
                 cuDoubleComplex *da, magma_int_t ldda,
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

    as returned by ZGEQRF_GPU.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix Q. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix Q. M >= N >= 0.

    K       (input) INTEGER
            The number of elementary reflectors whose product defines the
            matrix Q. N >= K >= 0.

    DA      (input/output) COMPLEX_16 array A on the GPU device, 
            dimension (LDDA,N). On entry, the i-th column must contain 
            the vector which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by ZGEQRF_GPU in the first k 
            columns of its array argument A.
            On exit, the M-by-N matrix Q.

    LDDA    (input) INTEGER
            The first dimension of the array A. LDDA >= max(1,M).

    TAU     (input) COMPLEX_16 array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZGEQRF_GPU.

    DT      (input) COMPLEX_16 work space array on the GPU device.
            This must be the 6th argument of magma_zgeqrf_gpu.

    NB      (input) INTEGER
            This is the block size used in ZGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and 
            stored in DT. 

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument has an illegal value
    =====================================================================    */

    #define da_ref(a_1,a_2) (da+(a_2)*(ldda) + (a_1))
    #define t_ref(a_1)      (dT+(a_1)*nb)

    magma_int_t  i__1, i__2, i__3;
    magma_int_t lwork;
    static magma_int_t i, ib, ki, kk, iinfo;
    magma_int_t lddwork = min(m, n);
    cuDoubleComplex *work, *panel;
    static cudaStream_t stream[2];

    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if ((n < 0) || (n > m)) {
	*info = -2;
    } else if ((k < 0) || (k > n)) {
	*info = -3;
    } else if (ldda < max(1,m)) {
	*info = -5;
    }
    if (*info != 0)
      return MAGMA_ERR_ILLEGAL_VALUE;

    if (n <= 0)
      return MAGMA_SUCCESS;

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

    /* Allocate work space on CPU in pinned memory */
    lwork = (n+m) * nb;
    if (kk < n)
      lwork = max(lwork, n * nb + (m-kk)*(n-kk));

    if ( cudaSuccess != 
         cudaMallocHost( (void**)&work, (lwork)*sizeof(cuDoubleComplex) ) )
      {
	*info = -11;
	return MAGMA_ERR_HOSTALLOC;
      }
    panel = work + n * nb;

    cudaStreamCreate(&stream[0]);
    cudaStreamCreate(&stream[1]);

    /* Use unblocked code for the last or only block. */
    if (kk < n)
      {
	i__1 = m - kk;
	i__2 = n - kk;
	i__3 = k - kk;
	cublasGetMatrix(i__1, i__2, sizeof(cuDoubleComplex),
			da_ref(kk, kk), ldda, panel, i__1);
	lapackf77_zungqr(&i__1, &i__2, &i__3, panel, &i__1, &tau[kk], 
			 work, &lwork, &iinfo);

	cublasSetMatrix(i__1, i__2, sizeof(cuDoubleComplex),
			panel, i__1, da_ref(kk, kk), ldda);
      }

    if (kk > 0)
      {
	/* Use blocked code */
	for (i = ki; i >= 0; i-=nb)
	  {
	    ib = min(nb, k - i);

	    /* Send current panel to the CPU for update */
	    i__2 = m - i;
	    cudaMemcpy2DAsync(panel,       i__2 * sizeof(cuDoubleComplex),
			      da_ref(i,i), ldda * sizeof(cuDoubleComplex),
			      sizeof(cuDoubleComplex)*i__2, ib,
			      cudaMemcpyDeviceToHost,stream[0]);

	    if (i + ib < n)
	      {
		/* Apply H to A(i:m,i+ib:n) from the left */
		i__3 = n - i - ib;
		magma_zlarfb_gpu( MagmaLeft, MagmaNoTrans, MagmaForward, MagmaColumnwise,
				  i__2, i__3, ib,
				  da_ref(i, i   ), ldda, t_ref(i),             nb,
				  da_ref(i, i+ib), ldda, dT + 2*lddwork*nb, lddwork);
	      }

	    /* Apply H to rows i:m of current block on the CPU */
	    cudaStreamSynchronize(stream[0]);
	    lapackf77_zungqr(&i__2, &ib, &ib, panel, &i__2, &tau[i], 
			     work, &lwork, &iinfo);
	    cudaMemcpy2DAsync(da_ref(i,i), ldda * sizeof(cuDoubleComplex),
			      panel,       i__2 * sizeof(cuDoubleComplex),
			      sizeof(cuDoubleComplex)*i__2, ib,
                              cudaMemcpyHostToDevice,stream[1]);

	    /* Set rows 1:i-1 of current block to zero */
            i__2 = i + ib;
	    magmablas_zlaset(i, i__2 - i, da_ref(0,i), ldda);
	  }
      }

    cudaFreeHost(work);
    cudaStreamDestroy(stream[0]);
    cudaStreamDestroy(stream[1]);

    return MAGMA_SUCCESS;
} /* magma_zungqr_gpu */

#undef da_ref
#undef t_ref
