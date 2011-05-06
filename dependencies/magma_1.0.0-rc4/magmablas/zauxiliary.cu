/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/
#include "common_magma.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- This is an auxiliary routine called from zgehrd.  The routine is called
      in 16 blocks, 32 thread per block and initializes to zero the 1st 
      32x32 block of A.
*/

__global__ void zset_to_zero(cuDoubleComplex *A, int lda){
    int ind = blockIdx.x*lda + threadIdx.x;
    
    A += ind;
    A[0] = MAGMA_Z_ZERO;
//   A[16*lda] = 0.;
}

__global__ void zset_nbxnb_to_zero(int nb, cuDoubleComplex *A, int lda){
   int ind = blockIdx.x*lda + threadIdx.x, i, j;

   A += ind;
   for(i=0; i<nb; i+=32){
     for(j=0; j<nb; j+=32)
         A[j] = MAGMA_Z_ZERO;
     A += 32*lda;
   }
}

void zzero_32x32_block(cuDoubleComplex *A, int lda)
{
  // zset_to_zero<<<16, 32>>>(A, lda);
  zset_to_zero<<<32, 32>>>(A, lda);
}

void zzero_nbxnb_block(int nb, cuDoubleComplex *A, int lda)
{
  zset_nbxnb_to_zero<<<32, 32>>>(nb, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- GPU kernel for initializing a matrix by 0
*/
#define zlaset_threads 64

__global__ void zlaset(int m, int n, cuDoubleComplex *A, int lda){
   int ibx = blockIdx.x * zlaset_threads;
   int iby = blockIdx.y * 32;

   int ind = ibx + threadIdx.x;

   A += ind + __mul24(iby, lda);

   #pragma unroll
   for(int i=0; i<32; i++)
     if (iby+i < n && ind < m)
        A[i*lda] = MAGMA_Z_ZERO;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Set the m x n matrix pointed by A to 0 on the GPU.
*/
extern "C" void
magmablas_zlaset(magma_int_t m, magma_int_t n, 
                 cuDoubleComplex *A, magma_int_t lda)
{
   dim3 threads(zlaset_threads, 1, 1);
   dim3 grid(m/zlaset_threads+(m % zlaset_threads != 0), n/32+(n%32!=0));

   zlaset<<< grid, threads >>> (m, n, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Given two matrices, 'a' on the CPU and 'da' on the GPU, this function
      returns the Frobenious norm of the difference of the two matrices.
      The function is used for debugging.
*/
double cpu_gpu_zdiff(int M, int N, cuDoubleComplex * a, int lda, cuDoubleComplex *da, int ldda)
{
  int one = 1, j;
  cuDoubleComplex mone = MAGMA_Z_NEG_ONE;
  double  work[1];
  cuDoubleComplex *ha = (cuDoubleComplex*)malloc( M * N * sizeof(cuDoubleComplex));
  double res;

  cublasGetMatrix(M, N, sizeof(cuDoubleComplex), da, ldda, ha, M);
  for(j=0; j<N; j++)
    blasf77_zaxpy(&M, &mone, a+j*lda, &one, ha+j*M, &one);
  res = lapackf77_zlange("f", &M, &N, ha, &M, work);

  free(ha);
  return res;
}

