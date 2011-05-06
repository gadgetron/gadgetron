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
   -- This is an auxiliary routine called from sgehrd.  The routine is called
      in 16 blocks, 32 thread per block and initializes to zero the 1st 
      32x32 block of A.
*/

__global__ void sset_to_zero(float *A, int lda){
    int ind = blockIdx.x*lda + threadIdx.x;
    
    A += ind;
    A[0] = MAGMA_S_ZERO;
//   A[16*lda] = 0.;
}

__global__ void sset_nbxnb_to_zero(int nb, float *A, int lda){
   int ind = blockIdx.x*lda + threadIdx.x, i, j;

   A += ind;
   for(i=0; i<nb; i+=32){
     for(j=0; j<nb; j+=32)
         A[j] = MAGMA_S_ZERO;
     A += 32*lda;
   }
}

void szero_32x32_block(float *A, int lda)
{
  // sset_to_zero<<<16, 32>>>(A, lda);
  sset_to_zero<<<32, 32>>>(A, lda);
}

void szero_nbxnb_block(int nb, float *A, int lda)
{
  sset_nbxnb_to_zero<<<32, 32>>>(nb, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- GPU kernel for initializing a matrix by 0
*/
#define slaset_threads 64

__global__ void slaset(int m, int n, float *A, int lda){
   int ibx = blockIdx.x * slaset_threads;
   int iby = blockIdx.y * 32;

   int ind = ibx + threadIdx.x;

   A += ind + __mul24(iby, lda);

   #pragma unroll
   for(int i=0; i<32; i++)
     if (iby+i < n && ind < m)
        A[i*lda] = MAGMA_S_ZERO;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Set the m x n matrix pointed by A to 0 on the GPU.
*/
extern "C" void
magmablas_slaset(magma_int_t m, magma_int_t n, 
                 float *A, magma_int_t lda)
{
   dim3 threads(slaset_threads, 1, 1);
   dim3 grid(m/slaset_threads+(m % slaset_threads != 0), n/32+(n%32!=0));

   slaset<<< grid, threads >>> (m, n, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Given two matrices, 'a' on the CPU and 'da' on the GPU, this function
      returns the Frobenious norm of the difference of the two matrices.
      The function is used for debugging.
*/
float cpu_gpu_sdiff(int M, int N, float * a, int lda, float *da, int ldda)
{
  int one = 1, j;
  float mone = MAGMA_S_NEG_ONE;
  float  work[1];
  float *ha = (float*)malloc( M * N * sizeof(float));
  float res;

  cublasGetMatrix(M, N, sizeof(float), da, ldda, ha, M);
  for(j=0; j<N; j++)
    blasf77_saxpy(&M, &mone, a+j*lda, &one, ha+j*M, &one);
  res = lapackf77_slange("f", &M, &N, ha, &M, work);

  free(ha);
  return res;
}

