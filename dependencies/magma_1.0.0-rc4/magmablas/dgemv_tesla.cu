/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"

#define magmablas_dgemv_tesla magmablas_dgemv
#define magmablas_dgemvt_tesla magmablas_dgemvt

#define num_threads 64
#define dgemv_bs 64

__global__ void 
dgemv_kernel(magma_int_t n, magma_int_t m, magma_int_t n1, double* A, magma_int_t lda, double *x, double *y)
{
  magma_int_t ind = blockIdx.x*num_threads + threadIdx.x;

  A += ind;
  x += threadIdx.x;

  double res = 0.f;

  __shared__ double buff[dgemv_bs];
  for(magma_int_t i=0; i<n1; i += dgemv_bs ){
    __syncthreads();
    buff[threadIdx.x]  = x[i];

    __syncthreads();
    #pragma unroll
    for(magma_int_t j=0; j < dgemv_bs ; j++){
       res+=A[0]*buff[j];
       A+=lda;
    }
  }
  __syncthreads();

  if (m>n1){
     buff[threadIdx.x]  = x[n1];

     __syncthreads();
     for(magma_int_t j=0; j<(m-n1); j++){
         res += A[0]*buff[j];
         A+=lda;
     }
  }

  if (ind<n)
     y[ind] = res;
}

extern "C" void
magmablas_dgemvt_tesla(magma_int_t m, magma_int_t n, double alpha, double *A, magma_int_t lda,
                       double *x, double *z);

extern "C" void
magmablas_dgemv_tesla(char trans,
                      magma_int_t m, magma_int_t n,
                      double alpha, 
                      double *A, magma_int_t lda, 
                      double *x, magma_int_t incx,
                      double beta,
                      double *z, magma_int_t incz)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    This routine computes:
    1) z =       A   x    if trans == 'N' or 'n', alpha == 1, beta == 0,
                          and incx == incz == 1 (using magmablas code)
    2) z = alpha A^t x    if trans == 'T' or 't', beta == 0,
                          and incx == incz == 1 (using magmablas code)
    3) z = alpha A^trans x + beta z
                          otherwise, using CUBLAS.

   Arguments
   ==========
    TRANS  - CHARACTER*1
             On entry, TRANS specifies the operation to be performed as
             follows:
               TRANS = 'N' or 'n'   z := alpha*A *x + beta*z
               TRANS = 'T' or 't'   z := alpha*A'*x + beta*z

    M      - (input) INTEGER
             On entry, M specifies the number of rows of the matrix A.

    N      - (input) INTEGER.
             On entry, N specifies the number of columns of the matrix A

    ALPHA  - DOUBLE PRECISION
             On entry, ALPHA specifies the scalar alpha.
             Unchanged on exit.

    A      - (input) DOUBLE PRECISION array of dimension ( LDA, n ) on the GPU.
   
    LDA    - (input) INTEGER
             LDA specifies the leading dimension of A.

    X      - (input) DOUBLE PRECISION array of dimension 
             n if trans == 'n'
             m if trans == 't'
      
    INCX   - (input) Specifies the increment for the elements of X.
             INCX must not be zero. Unchanged on exit.

    BETA   - DOUBLE PRECISION
             On entry, BETA specifies the scalar beta. When BETA is
             supplied as zero then Y need not be set on input.
             Unchanged on exit

    Z      - (output) DOUBLE PRECISION array of	dimension 
             m if trans == 'n'
             n if trans == 't' 
             
    INCZ  - (input) Specifies the increment for the elements of Z.
            INCZ must not be zero. Unchanged on exit.
    ===================================================================== */

    if (incx == 1 && incz == 1 && beta == 0.)
       if (trans == 'n' || trans == 'N')
          if (alpha == 1.)
            {
              magma_int_t blocks;
              if (m % num_threads==0)
                 blocks = m/num_threads;
              else
                 blocks = m/num_threads + 1;

              dim3 grid(blocks, 1, 1);
              dim3 threads(num_threads, 1, 1);
 
              dgemv_kernel<<<grid, threads>>>(m, n, 
                                              (n/dgemv_bs)*dgemv_bs, 
                                              A, lda, x, z);
          }
         else
           cublasDgemv(trans, m, n, alpha, A, lda,
                       x, incx, beta, z, incz);
       else
         magmablas_dgemvt_tesla(m, n, alpha, A, lda, x, z);
    else
      cublasDgemv(trans, m, n, alpha, A, lda,
                  x, incx, beta, z, incz);

}

__global__ void
dgemvt_kernel1(magma_int_t n, magma_int_t m, double alpha, magma_int_t n1, double* A, magma_int_t lda,
              double *x, double *y)
{
  const magma_int_t inx = threadIdx.x;
  const magma_int_t iny = threadIdx.y;

  magma_int_t ind  = iny + __mul24(blockIdx.x,32);
  ind = inx + __mul24(ind,lda);
  magma_int_t ind2 = inx + __mul24(iny,32);

  A += ind;
  x += ind2;

  double res = 0.f;

  __shared__ double buff[dgemv_bs];
  __shared__ double la[32][33];

  for(magma_int_t i=0; i<n1; i += dgemv_bs ){
      buff[ind2]  = x[i];
      #pragma unroll
      for(magma_int_t j=0; j<16; j++)
         la[iny+__mul24(j,2)][inx] = A[j*__mul24(2,lda)];

      __syncthreads();
      #pragma unroll
      for(magma_int_t j=0; j < 16; j++)
        res += la[inx][iny*16+j]*buff[j+iny*16];

      A += 32;

      //===============================================
      #pragma unroll
      for(magma_int_t j=0; j<16; j++)
         la[iny+__mul24(j,2)][inx] = A[j*__mul24(2,lda)];

      __syncthreads();

      #pragma unroll
      for(magma_int_t j=0; j < 16; j++)
        res += la[inx][iny*16+j]*buff[j+32+iny*16];
      A += 32;
    }

    if (n>n1){
      if (ind2>=(n-n1))
         buff[ind2]=0.;
      else
         buff[ind2]  = x[n1];

      #pragma unroll
      for(magma_int_t j=0; j<16; j++)
         la[iny+__mul24(j,2)][inx] = A[j*__mul24(2,lda)];

     __syncthreads();

     if (n-n1>16){
	#pragma unroll
        for(magma_int_t j=0; j < 16; j++)
           res += la[inx][iny*16+j]*buff[j+iny*16];

        A += 32;
        #pragma unroll
        for(magma_int_t j=0; j<16; j++)
          la[iny+__mul24(j,2)][inx] = A[j*__mul24(2,lda)];

	__syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 16; j++)
           res += la[inx][iny*16+j]*buff[j+32+iny*16];
     }
     else {
        #pragma unroll
        for(magma_int_t j=0; j < 16; j++)
          res += la[inx][iny*16+j]*buff[j+iny*16];
     }
  }
  ind = inx + __mul24(blockIdx.x,32);

  la[inx][iny]= res;
  if (ind<m){
     res = la[inx][0] + la[inx][1];
     y[ind] = alpha*res;
  }
}

__global__ void
dgemvt_kernel2(magma_int_t n, magma_int_t m, double alpha,
		magma_int_t n1, double* A, magma_int_t lda, double *x, double *y)
{
  const magma_int_t inx = threadIdx.x;
  const magma_int_t iny = threadIdx.y;

  magma_int_t ind  = iny + __mul24(blockIdx.x,16);
  ind = inx + __mul24(ind,lda);
  magma_int_t ind2 = inx + __mul24(iny,16);
  if (ind2>31)
     ind2-=32;

  A += ind;
  x += ind2;
  if (ind2>31)
     ind2-=32;

  double res = 0.f;

  __shared__ double buff[32];
  __shared__ double la[16][17];

  for(magma_int_t i=0; i<n1; i += 32 ){
     buff[ind2]  = x[i];
     #pragma unroll
     for(magma_int_t j=0; j<4; j++)
        la[iny+__mul24(j,4)][inx] = A[j*__mul24(4,lda)];

     __syncthreads();
     #pragma unroll
     for(magma_int_t j=0; j < 4; j++)
       res += la[inx][iny*4+j]*buff[j+iny*4];

     A += 16;
    __syncthreads();
     //===========================================
     #pragma unroll
     for(magma_int_t j=0; j<4; j++)
         la[iny+__mul24(j,4)][inx] = A[j*__mul24(4,lda)];

     __syncthreads();

     #pragma unroll
     for(magma_int_t j=0; j < 4; j++)
        res += la[inx][iny*4+j]*buff[j+16+iny*4];
     A += 16;
  }

  if (n>n1){
     if (ind2>=(n-n1))
        buff[ind2]=0.;
     else
        buff[ind2]  = x[n1];

     __syncthreads();
     #pragma unroll
     for(magma_int_t j=0; j<4; j++)
         la[iny+__mul24(j,4)][inx] = A[j*__mul24(4,lda)];

     __syncthreads();
     if (n-n1>4){
        #pragma unroll
	for(magma_int_t j=0; j < 4; j++)
           res += la[inx][iny*4+j]*buff[j+iny*4];

        A += 16;
        __syncthreads();
        #pragma unroll
          for(magma_int_t j=0; j<4; j++)
            la[iny+__mul24(j,4)][inx] = A[j*__mul24(4,lda)];

        __syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 4; j++)
           res += la[inx][iny*4+j]*buff[j+16+iny*4];
     }
     else {
        #pragma unroll
        for(magma_int_t j=0; j < 4; j++)
          res += la[inx][iny*4+j]*buff[j+iny*4];
     }
  }

  __syncthreads();
  ind = inx + __mul24(blockIdx.x,16);
  la[inx][iny]= res;
  __syncthreads();
  if (ind<m){
     res = la[inx][0] + la[inx][1] + la[inx][2] + la[inx][3];
     y[ind] = alpha*res;
  }
}

extern "C" void
magmablas_dgemvt1_tesla(magma_int_t m, magma_int_t n, double alpha, double *A, magma_int_t lda,
                        double *x, double *z)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    This routine computes z = alpha A^t x on the GPU.
    Recommended for large M and N.

    M      - (input) INTEGER.
             On entry, M specifies the number of rows of the matrix A.

    N      - (input) INTEGER.
             On entry, N specifies the number of columns of the matrix A

    A      - (input) DOUBLE PRECISION array of dimension ( LDA, n ) on the GPU.

    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) DOUBLE PRECISION array of dimension m.

    Z      - (output) DOUBLE PRECISION array of dimension n.
             On exit Z = alpha A^t X.
    ===================================================================== */

    magma_int_t blocks;

    if (n % 32==0)
        blocks = n/32;
    else
        blocks = n/32 + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(32, 2, 1);
    dgemvt_kernel1<<<grid, threads>>>(m, n, alpha, (m / dgemv_bs)*dgemv_bs,
                                      A, lda, x, z);
}

extern "C" void
magmablas_dgemvt2_tesla(magma_int_t m, magma_int_t n, double alpha, double *A, magma_int_t lda,
                        double *x, double *z)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    This routine computes z = alpha A^t x on the GPU. Used in least squares
    solver for N small (e.g. = BS, a block size of order 64, 128, etc).

    M      - (input) INTEGER.
             On entry, M specifies the number of rows of the matrix A.

    N      - (input) INTEGER.
             On entry, N specifies the number of columns of the matrix A

    A      - (input) DOUBLE PRECISION array of dimension ( LDA, n ) on the GPU.

    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) DOUBLE PRECISION array of dimension m.

    Z      - (output) DOUBLE PRECISION array of dimension n.
             On exit Z = alpha A^t X.
    ===================================================================== */

    magma_int_t blocks;

    if (n % 16==0)
        blocks = n/16;
    else
        blocks = n/16 + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(16, 4, 1);

    dgemvt_kernel2<<<grid, threads>>>(m, n, alpha, (m / 32)*32,
                                      A, lda, x, z);
}

extern "C" void
magmablas_dgemvt_tesla(magma_int_t m, magma_int_t n, double alpha, double *A, magma_int_t lda,
                       double *x, double *z)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    This routine computes z = alpha A^t x on the GPU.

    M      - (input) INTEGER.
             On entry, M specifies the number of rows of the matrix A.

    N      - (input) INTEGER.
             On entry, N specifies the number of columns of the matrix A

    A      - (input) SINGLE PRECISION array of dimension (LDA, n) on the GPU.

    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) SINGLE PRECISION array of dimension m.

    Z      - (output) SINGLE PRECISION array of dimension n.
             On exit Z = alpha A^t X.
    ===================================================================== */

    if (n<=128)
      magmablas_dgemvt2_tesla(m, n, alpha, A, lda, x, z);
    else
      magmablas_dgemvt1_tesla(m, n, alpha, A, lda, x, z);
}

#undef num_threads
#undef dgemv_bs
