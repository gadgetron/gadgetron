/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated s

*/
#include "common_magma.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Auxiliary function: 'a' is pointer to the current panel holding the 
      Householder vectors for the QR factorization of the panel. This routine
      puts ones on the diagonal and zeros in the upper triangular part of 'a'.
      The upper triangular values are stored in work. Than the inverse is 
      calculated in place in work, so as final result work holds the inverse
      of the upper triangular diagonal block.
 */
void ssplit_diag_block(int ib, float *a, int lda, float *work){
    int i, j, info;
    float *cola, *colw;
    float c_zero = MAGMA_S_ZERO;
    float c_one  = MAGMA_S_ONE;

    for(i=0; i<ib; i++){
        cola = a    + i*lda;
        colw = work + i*ib;
        for(j=0; j<i; j++){
            colw[j] = cola[j];
            cola[j] = c_zero;
        }
        colw[i] = cola[i];
        cola[i] = c_one;
    }
    lapackf77_strtri( MagmaUpperStr, MagmaNonUnitStr, &ib, work, &ib, &info);
}

extern "C" magma_int_t
magma_sgeqrf_gpu( magma_int_t m, magma_int_t n, 
		  float *dA,   magma_int_t ldda,
                  float *tau, float *dT, 
		  magma_int_t *info )
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    SGEQRF computes a QR factorization of a REAL M-by-N matrix A:
    A = Q * R. This version stores the triangular matrices used in
    the factorization so that they can be applied directly (i.e.,
    without being recomputed) later. As a result, the application
    of Q is much faster.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) REAL array on the GPU, dimension (LDDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).

    LDDA     (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).
            To benefit from coalescent memory accesses LDDA must be
            dividable by 16.

    TAU     (output) REAL array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    dT      (workspace/output)  REAL array on the GPU, 
            dimension (2*MIN(M, N) + (N+31)/32*32 )*NB,
            where NB can be obtained through magma_get_sgeqrf_nb(M).
            It starts with MIN(M,N)*NB block that store the triangular T
            matrices, followed by the MIN(M,N)*NB block of the diagonal
            inverses for the R matrix. The rest of the array is used as workspace.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -9, internal GPU memory allocation failed.

    Further Details
    ===============
    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).
    =====================================================================    */

    #define a_ref(a_1,a_2) (dA+(a_2)*(ldda) + (a_1))
    #define t_ref(a_1)     (dT+(a_1)*nb)
    #define d_ref(a_1)     (dT+(lddwork+(a_1))*nb)
    #define dd_ref(a_1)    (dT+(2*lddwork+(a_1))*nb)
    #define work_ref(a_1)  ( work + (a_1))
    #define hwork          ( work + (nb)*(m))

    magma_int_t i, k, old_i, old_ib, rows, cols;
    magma_int_t ib, nb;
    magma_int_t ldwork, lddwork, lwork, lhwork;
    float *work, *ut;

    /* check arguments */
    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0)
        return MAGMA_ERR_ILLEGAL_VALUE;

    k = min(m,n);
    if (k == 0)
        return MAGMA_SUCCESS;

    nb = magma_get_sgeqrf_nb(m);

    lwork  = (m + n + nb)*nb;
    lhwork = lwork - m*nb;

    if ( cudaSuccess != cudaMallocHost((void**)&work, lwork*sizeof(float)) ) {
	*info = -9;
	magma_xerbla( "magma_sgeqrf_gpu", info );
	return MAGMA_ERR_HOSTALLOC;
    }
    
    ut = hwork+nb*(n);
    memset( ut, 0, nb*nb*sizeof(float));

    static cudaStream_t stream[2];
    cudaStreamCreate(&stream[0]);
    cudaStreamCreate(&stream[1]);

    ldwork = m;
    lddwork= k;

    if ( (nb > 1) && (nb < k) ) {
        /* Use blocked code initially */
        old_i = 0; old_ib = nb;
        for (i = 0; i < k-nb; i += nb) {
            ib = min(k-i, nb);
            rows = m -i;
            cudaMemcpy2DAsync( work_ref(i), ldwork*sizeof(float),
                               a_ref(i,i),  ldda   *sizeof(float),
                               sizeof(float)*rows, ib,
                               cudaMemcpyDeviceToHost, stream[1]);
            if (i>0){
                /* Apply H' to A(i:m,i+2*ib:n) from the left */
                cols = n-old_i-2*old_ib;
                magma_slarfb_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
				  m-old_i, cols, old_ib,
				  a_ref(old_i, old_i         ), ldda, t_ref(old_i), nb,
				  a_ref(old_i, old_i+2*old_ib), ldda, dd_ref(0),    lddwork);
		
                /* store the diagonal */
                cudaMemcpy2DAsync(d_ref(old_i), old_ib * sizeof(float),
                                  ut,           old_ib * sizeof(float),
                                  sizeof(float)*old_ib, old_ib,
                                  cudaMemcpyHostToDevice, stream[0]);
            }

            cudaStreamSynchronize(stream[1]);
            lapackf77_sgeqrf(&rows, &ib, work_ref(i), &ldwork, tau+i, hwork, &lhwork, info);
            /* Form the triangular factor of the block reflector
               H = H(i) H(i+1) . . . H(i+ib-1) */
            lapackf77_slarft( MagmaForwardStr, MagmaColumnwiseStr, 
                              &rows, &ib, 
                              work_ref(i), &ldwork, tau+i, hwork, &ib);

            /* Put 0s in the upper triangular part of a panel (and 1s on the
               diagonal); copy the upper triangular in ut and invert it     */
            cudaStreamSynchronize(stream[0]);
            ssplit_diag_block(ib, work_ref(i), ldwork, ut);
            cublasSetMatrix(rows, ib, sizeof(float),
                            work_ref(i), ldwork, a_ref(i,i), ldda);

            if (i + ib < n) {
                /* Send the triangular factor T to the GPU */
                cublasSetMatrix(ib, ib, sizeof(float), hwork, ib, t_ref(i), nb);

                if (i+nb < k-nb){
                    /* Apply H' to A(i:m,i+ib:i+2*ib) from the left */
                    magma_slarfb_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
				      rows, ib, ib, 
				      a_ref(i, i   ), ldda, t_ref(i),  nb, 
				      a_ref(i, i+ib), ldda, dd_ref(0), lddwork);
                }
                else {
                    cols = n-i-ib;
                    magma_slarfb_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
				      rows, cols, ib, 
				      a_ref(i, i   ), ldda, t_ref(i),  nb, 
				      a_ref(i, i+ib), ldda, dd_ref(0), lddwork);
                    /* Fix the diagonal block */
                    cublasSetMatrix(ib, ib, sizeof(float), ut, ib, d_ref(i), ib);
                }
                old_i  = i;
                old_ib = ib;
            }
        }
    } else {
        i = 0;
    }

    /* Use unblocked code to factor the last or only block. */
    if (i < k) {
        ib   = n-i;
        rows = m-i;
        cublasGetMatrix(rows, ib, sizeof(float),
                        a_ref(i, i), ldda, 
                        work,        rows);
        lhwork = lwork - rows*ib;
        lapackf77_sgeqrf(&rows, &ib, work, &rows, tau+i, work+ib*rows, &lhwork, info);
        
        cublasSetMatrix(rows, ib, sizeof(float),
                        work,        rows, 
                        a_ref(i, i), ldda);
    }

    cudaStreamDestroy(stream[0]);
    cudaStreamDestroy(stream[1]);
    cudaFreeHost(work);
    return MAGMA_SUCCESS;

/*     End of MAGMA_SGEQRF */

} /* magma_sgeqrf */

#undef a_ref
#undef t_ref
#undef d_ref
#undef work_ref
