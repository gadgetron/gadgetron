/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_d
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define cublasDgemm magmablas_dgemm
  #define cublasDtrsm magmablas_dtrsm
#endif

#if (GPUSHMEM >= 200)
  #if (defined(PRECISION_s))
     #undef  cublasSgemm
     #define cublasSgemm magmablas_sgemm_fermi80
  #endif
#endif
// === End defining what BLAS to use =======================================

#define dA(i, j)  (dA + (j)*ldda + (i))

extern "C" magma_int_t
magma_dpotrf_gpu(char uplo, magma_int_t n, 
                 double *dA, magma_int_t ldda, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   
    DPOTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix dA.   

    The factorization has the form   
       dA = U\*\*H * U,  if UPLO = 'U', or   
       dA = L  * L\*\*H,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   
    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of dA is stored;   
            = 'L':  Lower triangle of dA is stored.   

    N       (input) INTEGER   
            The order of the matrix dA.  N >= 0.   

    dA      (input/output) DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)   
            On entry, the symmetric matrix dA.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of dA contains the upper   
            triangular part of the matrix dA, and the strictly lower   
            triangular part of dA is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of dA contains the lower   
            triangular part of the matrix dA, and the strictly upper   
            triangular part of dA is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization dA = U\*\*H*U or dA = L*L\*\*H.   

    LDDA     (input) INTEGER   
            The leading dimension of the array dA.  LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            dividable by 16.

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   
    =====================================================================   */


    magma_int_t     j, jb, nb;
    char            uplo_[2] = {uplo, 0};
    double zone  = MAGMA_D_ONE;
    double mzone = MAGMA_D_NEG_ONE;
    double *work;
    double          done  = (double) 1.0;
    double          mdone = (double)-1.0;
    long int        upper = lapackf77_lsame(uplo_, "U");

    *info = 0;
    if ( (! upper) && (! lapackf77_lsame(uplo_, "L")) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0)
        return MAGMA_ERR_ILLEGAL_VALUE;

    nb = magma_get_dpotrf_nb(n);

    if (cudaSuccess != cudaMallocHost( (void**)&work, nb*nb*sizeof(double) ) ) {
	*info = -6;
	return MAGMA_ERR_HOSTALLOC;
    }

    static cudaStream_t stream[2];
    cudaStreamCreate(&stream[0]);
    cudaStreamCreate(&stream[1]);

    if ((nb <= 1) || (nb >= n)) {
        /*  Use unblocked code. */
        cublasGetMatrix(n, n, sizeof(double), dA, ldda, work, n);
        lapackf77_dpotrf(uplo_, &n, work, &n, info);
        cublasSetMatrix(n, n, sizeof(double), work, n, dA, ldda);
    } else {

        /* Use blocked code. */
	if (upper) {
            
            /* Compute the Cholesky factorization A = U'*U. */
            for (j=0; j<n; j+=nb) {
                
                /* Update and factorize the current diagonal block and test   
                   for non-positive-definiteness. Computing MIN */
		jb = min(nb, (n-j));
                
                cublasDsyrk(MagmaUpper, MagmaTrans, jb, j, 
                            mdone, dA(0, j), ldda, 
                            done,  dA(j, j), ldda);

                cudaMemcpy2DAsync(work,     jb  *sizeof(double), 
                                  dA(j, j), ldda*sizeof(double), 
                                  jb*sizeof(double), jb, 
                                  cudaMemcpyDeviceToHost,stream[1]);
		
		if ( (j+jb) < n) {
                    /* Compute the current block row. */
                    cublasDgemm(MagmaTrans, MagmaNoTrans, 
                                jb, (n-j-jb), j,
                                mzone, dA(0, j   ), ldda, 
                                       dA(0, j+jb), ldda,
                                zone,  dA(j, j+jb), ldda);
                }
                
                cudaStreamSynchronize(stream[1]);

                lapackf77_dpotrf(MagmaUpperStr, &jb, work, &jb, info);
                cudaMemcpy2DAsync( dA(j, j), ldda*sizeof(double), 
                                   work,     jb  *sizeof(double), 
                                   sizeof(double)*jb, jb, 
                                   cudaMemcpyHostToDevice,stream[0]);
		if (*info != 0) {
		  *info = *info + j;
		  break;
                }

                if ( (j+jb) < n)
                    cublasDtrsm( MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit, 
                                 jb, (n-j-jb),
                                 zone, dA(j, j   ), ldda, 
                                       dA(j, j+jb), ldda);
	    }
	} else {
            //=========================================================
            // Compute the Cholesky factorization A = L*L'.
	    for (j=0; j<n; j+=nb) {

                //  Update and factorize the current diagonal block and test   
                //  for non-positive-definiteness. Computing MIN 
                jb = min(nb, (n-j));

                cublasDsyrk(MagmaLower, MagmaNoTrans, jb, j,
                            mdone, dA(j, 0), ldda, 
                            done,  dA(j, j), ldda);
		
                cudaMemcpy2DAsync( work,     jb  *sizeof(double),
                                   dA(j, j), ldda*sizeof(double),
                                   sizeof(double)*jb, jb,
                                   cudaMemcpyDeviceToHost,stream[1]);
		
                if ( (j+jb) < n) {
                    cublasDgemm( MagmaNoTrans, MagmaTrans, 
                                 (n-j-jb), jb, j,
                                 mzone, dA(j+jb, 0), ldda, 
                                        dA(j,    0), ldda,
                                 zone,  dA(j+jb, j), ldda);
                }

                cudaStreamSynchronize(stream[1]);
	        lapackf77_dpotrf(MagmaLowerStr, &jb, work, &jb, info);
	        cudaMemcpy2DAsync(dA(j, j), ldda*sizeof(double), 
                                  work,     jb  *sizeof(double), 
                                  sizeof(double)*jb, jb, 
                                  cudaMemcpyHostToDevice,stream[0]);
		if (*info != 0) {
		  *info = *info + j;
		  break;
                }
	        
		if ( (j+jb) < n)
                    cublasDtrsm(MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit, 
                                (n-j-jb), jb, 
                                zone, dA(j,    j), ldda, 
                                      dA(j+jb, j), ldda);
	    }

	}
    }

    cudaStreamDestroy(stream[0]);
    cudaStreamDestroy(stream[1]);
    cudaFreeHost(work);

    return MAGMA_SUCCESS;
} /* magma_dpotrf_gpu */
