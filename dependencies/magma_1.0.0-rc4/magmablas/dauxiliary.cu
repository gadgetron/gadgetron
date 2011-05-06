/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d

*/
#include "common_magma.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- This is an auxiliary routine called from dgehrd.  The routine is called
      in 16 blocks, 32 thread per block and initializes to zero the 1st 
      32x32 block of A.
*/

__global__ void dset_to_zero(double *A, int lda){
    int ind = blockIdx.x*lda + threadIdx.x;
    
    A += ind;
    A[0] = MAGMA_D_ZERO;
//   A[16*lda] = 0.;
}

__global__ void dset_nbxnb_to_zero(int nb, double *A, int lda){
   int ind = blockIdx.x*lda + threadIdx.x, i, j;

   A += ind;
   for(i=0; i<nb; i+=32){
     for(j=0; j<nb; j+=32)
         A[j] = MAGMA_D_ZERO;
     A += 32*lda;
   }
}

void dzero_32x32_block(double *A, int lda)
{
  // dset_to_zero<<<16, 32>>>(A, lda);
  dset_to_zero<<<32, 32>>>(A, lda);
}

void dzero_nbxnb_block(int nb, double *A, int lda)
{
  dset_nbxnb_to_zero<<<32, 32>>>(nb, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- GPU kernel for initializing a matrix by 0
*/
#define dlaset_threads 64

__global__ void dlaset(int m, int n, double *A, int lda){
   int ibx = blockIdx.x * dlaset_threads;
   int iby = blockIdx.y * 32;

   int ind = ibx + threadIdx.x;

   A += ind + __mul24(iby, lda);

   #pragma unroll
   for(int i=0; i<32; i++)
     if (iby+i < n && ind < m)
        A[i*lda] = MAGMA_D_ZERO;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Set the m x n matrix pointed by A to 0 on the GPU.
*/
extern "C" void
magmablas_dlaset(magma_int_t m, magma_int_t n, 
                 double *A, magma_int_t lda)
{
   dim3 threads(dlaset_threads, 1, 1);
   dim3 grid(m/dlaset_threads+(m % dlaset_threads != 0), n/32+(n%32!=0));

   dlaset<<< grid, threads >>> (m, n, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Given two matrices, 'a' on the CPU and 'da' on the GPU, this function
      returns the Frobenious norm of the difference of the two matrices.
      The function is used for debugging.
*/
double cpu_gpu_ddiff(int M, int N, double * a, int lda, double *da, int ldda)
{
  int one = 1, j;
  double mone = MAGMA_D_NEG_ONE;
  double  work[1];
  double *ha = (double*)malloc( M * N * sizeof(double));
  double res;

  cublasGetMatrix(M, N, sizeof(double), da, ldda, ha, M);
  for(j=0; j<N; j++)
    blasf77_daxpy(&M, &mone, a+j*lda, &one, ha+j*M, &one);
  res = lapackf77_dlange("f", &M, &N, ha, &M, work);

  free(ha);
  return res;
}

