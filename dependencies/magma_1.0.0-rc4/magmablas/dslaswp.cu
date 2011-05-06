/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated ds

*/
#include "common_magma.h"

#define num_threadds 64

__global__ void
dslaswp_kernel(int n, double *a, int lda, float *sa, int m, int *ipiv)
{
    int ind = blockIdx.x*num_threadds + threadIdx.x, newind;
    float res;
    
    if (ind < m) {
        sa   += ind;
        ipiv += ind;
        
        newind = ipiv[0];
        
        for(int i=0; i<n; i++) {
            res = MAGMA_S_MAKE( (float)(a[newind+i*lda]),
                                (float)(a[newind+i*lda]) );
            sa[i*lda] = res; 
        }
    }
}

__global__ void
dslaswp_inv_kernel(int n, double *a, int lda, float *sa, int m, int *ipiv)
{
    int ind = blockIdx.x*num_threadds + threadIdx.x, newind;
    double res;

    if (ind < m) {
        a   += ind;
        ipiv += ind;

        newind = ipiv[0];

        for(int i=0; i<n; i++) {
            res = MAGMA_D_MAKE( (double)(sa[newind+i*lda]),
                                (double)(sa[newind+i*lda]) );
            a[i*lda] = res;
        }
    }
}


extern "C" void
magmablas_dslaswp(int n, double *a, int lda, float *sa, int m, int *ipiv, int incx)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    Row i of A is casted to single precision in row ipiv[i] of SA, 0<=i<m.

    N      - (input) INTEGER.
             On entry, N specifies the number of columns of the matrix A.

    A      - (input) DOUBLE PRECISION array on the GPU, dimension (LDA,N)
             On entry, the matrix of column dimension N and row dimension M
             to which the row interchanges will be applied.

    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    SA     - (output) REAL array on the GPU, dimension (LDA,N)
             On exit, the casted to single precision and permuted matrix.
        
    M      - (input) The number of rows to be interchanged.

    IPIV   - (input) INTEGER array, dimension (M)
             The vector of pivot indices. Row i of A is casted to single 
             precision in row ipiv[i] of SA, 0<=i<m. 

    INCX   - (input) INTEGER
             If IPIV is negative, the pivots are applied in reverse 
             order, otherwise in straight-forward order.
    ===================================================================== */

    int blocks;
    if (m % num_threadds==0)
	blocks = m/num_threadds;
    else
        blocks = m/num_threadds + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(num_threadds, 1, 1);

    if (incx >=0)
      dslaswp_kernel<<<grid, threads>>>(n, a, lda, sa, m, ipiv);
    else
      dslaswp_inv_kernel<<<grid, threads>>>(n, a, lda, sa, m, ipiv);
}

#undef num_threadds
