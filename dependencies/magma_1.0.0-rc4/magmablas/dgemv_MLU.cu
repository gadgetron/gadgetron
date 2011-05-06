/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"

#define num_threads 64
#define dgemv_bs 64

__global__ void 
magma_dgemv_MLU(int n, int m, int n1, double* A, int lda, double *x, double *y)
{
  int ind = blockIdx.x*num_threads + threadIdx.x;

  A += ind;
  x += threadIdx.x;

  double res = 0.f;

  __shared__ double buff[dgemv_bs];

  for(int i=0; i<n1; i += dgemv_bs ){
    __syncthreads();
    buff[threadIdx.x]  = x[i];

    __syncthreads();
    #pragma unroll
    for(int j=0; j < dgemv_bs ; j++){
       res+=A[0]*buff[j];
       A+=lda;
    }
  }
  __syncthreads();
  
  if (m>n1)
 {
   if( (threadIdx.x + n1 ) >= m ) {
	x+= ( m - threadIdx.x -1 ) ; 
   } 	
   else{
	x+=n1;
   }	
   m = m -n1 ; 
/*
	 Note
	 ====
	 Stan ............
	 This is going to give segmentation fault or Error in GPU for illegal memory access. -- I am talking about x index 
         buff[threadIdx.x]  = x[n1];
*/
     buff[threadIdx.x]  = x[0];

     __syncthreads();
     for(int j=0; j<m; j++){
         res += A[0]*buff[j];
         A+=lda;
     }
  }

  if (ind<n)
     y[ind] -= res;
}

extern "C" void
magmablas_dgemv_MLU(int n, int m, double *A, int lda, double *x, double *z)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    This routine computes z = z - Ax on the GPU.

    N      - (input) INTEGER.
             On entry, N specifies the number of rows of the matrix A.

    M      - (input) INTEGER.
             On entry, M specifies the number of columns of the matrix A

    A      - (input) DOUBLE PRECISION array of dimension ( LDA, m ) on the GPU.
   
    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) DOUBLE PRECISION array of dimension m.
     
    Z      - (output) DOUBLE PRECISION array of	dimension m. 
             On exit Z = A - A X.

    ===================================================================== */

    int blocks;
    if (n % num_threads==0)
        blocks = n/num_threads;
    else
        blocks = n/num_threads + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(num_threads, 1, 1);
 
    magma_dgemv_MLU <<<grid, threads>>>(n, m, (m / dgemv_bs)*dgemv_bs, 
                                    A, lda, x, z);
}

#undef num_threads
#undef dgemv_bs
