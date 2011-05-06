/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"

#define magmablas_sgemv_fermi magmablas_sgemv

#define num_threads 128
#define sgemv_bs 32
#define threadSize 128


__global__ void 
sgemvn_kernel1_fermi(magma_int_t n, magma_int_t m, magma_int_t n1, float alpha, float* A, magma_int_t lda, float *x, float *y)
{
  magma_int_t ind = blockIdx.x*num_threads + threadIdx.x;

  A += ind;

  float res = 0.f;

  for(magma_int_t i=0; i<n1; i += sgemv_bs ){

    #pragma unroll
    for(magma_int_t j=0; j < sgemv_bs ; j++){
       res += A[0] * x[j];
       A   += lda;
    }
	x += sgemv_bs;
  }

  if (m>n1){

     for(magma_int_t j=0; j<(m-n1); j++){
         res += A[0] * x[j];
         A   += lda;
     }
  }

  if (ind<n)
     y[ind] = alpha * res;

}

__global__ void 
sgemvn_kernel2_fermi(magma_int_t n, magma_int_t m, magma_int_t n1, float alpha,  float* A, magma_int_t lda, float *x, float *y)
{
  magma_int_t ind = blockIdx.x*num_threads + threadIdx.x;

  A += ind;
  x += threadIdx.x;

  float res = 0.f;

  __shared__ float buff[num_threads];
  for(magma_int_t i=0; i<n1; i += num_threads ){
    __syncthreads();
    buff[threadIdx.x]  = x[i];

    __syncthreads();
    #pragma unroll
    for(magma_int_t j=0; j < num_threads ; j++){
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
     y[ind] = alpha * res;
}

extern "C" void
magmablas_sgemvn_fermi(magma_int_t n, magma_int_t m, float alpha, float *A, magma_int_t lda, float *x, float *y)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    This routine computes Y = alpha A x on the GPU.

    N      - (input) INTEGER.
             On entry, N specifies the number of rows of the matrix A.

    M      - (input) INTEGER.
             On entry, M specifies the number of columns of the matrix A

    A      - (input) SINGLE PRECISION array of dimension ( LDA, m ) on the GPU.
   
    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) SINGLE PRECISION array of dimension m.
     
    Y      - (output) SINGLE PRECISION array of	dimension m. 
             On exit Y = alpha A X.

    ===================================================================== */

    magma_int_t blocks;
    if (n % num_threads==0)
        blocks = n/num_threads;
    else
        blocks = n/num_threads + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(num_threads, 1, 1);
    if(n<=8500) 
		sgemvn_kernel1_fermi<<<grid, threads>>>(n, m, (m / sgemv_bs)*sgemv_bs, 
			                           alpha, A, lda, x, y);
    else
		sgemvn_kernel2_fermi<<<grid, threads>>>(n, m, (m / num_threads)*num_threads, 
			                           alpha, A, lda, x, y);
}



__global__ void 
sgemvt_kernel1_fermi(magma_int_t m, magma_int_t n, float alpha, magma_int_t n1, float* A, magma_int_t lda,
              float *x, float *y)
{
	magma_int_t tx = threadIdx.x;

	__shared__ float sdata[threadSize];
	
	volatile float *smem;

	float res;
	res = 0.0f;
     
	for(magma_int_t i=0; i<n1; i+= threadSize)
	{
		res += A[tx + i + lda * blockIdx.y] * x[tx + i];
	}

	
	if(m > n1)
	{
		if( tx + n1 <  m )
		{
			res  += A[tx + n1 + lda *blockIdx.y] * x[tx + n1];
		}
		else 
		{
			res  += 0.0f;
		}
	}	

    sdata[tx] = res;
	__syncthreads();
    
    /*
	if(tx < 128) 
	{
		sdata[tx] += sdata[tx + 128];
	}
    __syncthreads();
	*/

	if(tx < 64) 
	{
		sdata[tx] += sdata[tx + 64];
	}
    __syncthreads();

	if(tx < 32)
	{
		smem = sdata;
		smem[tx] += smem[tx + 32];
		smem[tx] += smem[tx + 16];
		smem[tx] += smem[tx +  8];
		smem[tx] += smem[tx +  4];
		smem[tx] += smem[tx +  2];
		smem[tx] += smem[tx +  1];
	}

    if( tx == 0 ) 
	{
		y[blockIdx.y] = sdata[0]; 		

		if (blockIdx.y < n)
		{
			y[blockIdx.y] = y[blockIdx.y] * alpha;
		}
	}
}


__global__ void 
sgemvt_kernel2_fermi(magma_int_t m, magma_int_t n, float alpha,
               magma_int_t n1, float* A, magma_int_t lda, float *x, float *y)
{
  const magma_int_t inx = threadIdx.x;
  const magma_int_t iny = threadIdx.y;

  magma_int_t ind  = iny + blockIdx.x * 16;
  ind = inx + ind * lda;
  magma_int_t ind2 = inx + iny * 16;
  if (ind2>31)
     ind2-=32;

  A += ind;
  x += ind2;

  float res = 0.f;

  __shared__ float buff[32];
  __shared__ float la[16][17];

  for(magma_int_t i=0; i<n1; i += 32 ){
     buff[ind2]  = x[i];
     #pragma unroll
     for(magma_int_t j=0; j<4; j++)
        la[iny + j * 4][inx] = A[j* 4 * lda];

     __syncthreads();
     #pragma unroll
     for(magma_int_t j=0; j < 4; j++)
       res += la[inx][iny*4+j]*buff[j+iny*4];

     A += 16;

     __syncthreads();
     //===========================================
     #pragma unroll
     for(magma_int_t j=0; j<4; j++)
         la[iny+ j * 4][inx] = A[j* 4 * lda];

     __syncthreads();

     #pragma unroll
     for(magma_int_t j=0; j < 4; j++)
        res += la[inx][iny*4+j]*buff[j+16+iny*4];
     A += 16;
  }

  __syncthreads(); // 1
  if (n>n1){
     if (ind2>=(n-n1))
	buff[ind2]=0.;
     else
        buff[ind2]  = x[n1];

     __syncthreads();
     #pragma unroll
     for(magma_int_t j=0; j<4; j++)
         if (inx>=(n-n1))
            la[iny + j * 4][inx] =  0.f;
         else
            la[iny + j * 4][inx] = A[j* 4 * lda];

     __syncthreads();
     if (n-n1>4){
        #pragma unroll
        for(magma_int_t j=0; j < 4; j++){
           ind =  j+iny*4;
           res += la[inx][ind]*buff[ind];
        }
	A += 16;
        __syncthreads();
	#pragma unroll
	for(magma_int_t j=0; j<4; j++)
          if (inx+16>=(n-n1))
             la[iny+ j * 4][inx] = 0.f;
          else
             la[iny+ j * 4][inx] = A[j* 4* lda];

        __syncthreads();

        #pragma unroll
	for(magma_int_t j=0; j < 4; j++){
           ind = j+4*iny;
           res += la[inx][ind]*buff[16+ind];
        }
     }
     else {
	#pragma unroll
        for(magma_int_t j=0; j < 4; j++){
          ind = j+iny*4;
          res += la[inx][ind]*buff[ind];
        }
     }
  }

  __syncthreads();
  ind = inx + blockIdx.x * 16;
  la[inx][iny]= res;
  __syncthreads();
  if (ind<n && iny==0){
     res = la[inx][0] + la[inx][1] + la[inx][2] + la[inx][3];
     y[ind] = alpha*res;
  }
}

extern "C" void
magmablas_sgemvt1_fermi(magma_int_t m, magma_int_t n, float alpha, float *A, magma_int_t lda,
                  float *x, float *y)
{


    dim3 grid    ( 1,  n,  1);
    dim3 threads ( threadSize,   1,  1);

    sgemvt_kernel1_fermi<<<grid, threads>>>( m, n, alpha, ( m / threadSize) * threadSize,
                                       A, lda, x, y);

									  
}

extern "C" void
magmablas_sgemvt2_fermi(magma_int_t m, magma_int_t n, float alpha, float *A, magma_int_t lda,
                  float *x, float *y)
{

    magma_int_t blocks;

    if (n % 16==0)
        blocks = n/16;
    else
        blocks = n/16 + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(16, 4, 1);

    sgemvt_kernel2_fermi<<<grid, threads>>>(m, n, alpha, (m / 32)*32,
                                      A, lda, x, y);
}

extern "C" void
magmablas_sgemvt_fermi(magma_int_t m, magma_int_t n, float alpha, float *A, magma_int_t lda, 
                 float *x, float *y)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    This routine computes y = alpha *  A^t *  x on the GPU.

    M      - (input) INTEGER.
             On entry, M specifies the number of rows of the matrix A.

    N      - (input) INTEGER.
             On entry, N specifies the number of columns of the matrix A

    A      - (input) SINGLE PRECISION array of dimension ( LDA, n ) on the GPU.

    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) SINGLE PRECISION array of dimension m.

    Y      - (output) SINGLE PRECISION array of dimension n.
             On exit Y = alpha A^t X.

    ===================================================================== */

    if (n<=128)
      magmablas_sgemvt2_fermi(m, n, alpha, A, lda, x, y);
    else
      magmablas_sgemvt1_fermi(m, n, alpha, A, lda, x, y);
    

}


extern "C" void
magmablas_sgemv_fermi(char trans,
                      magma_int_t m, magma_int_t n,
                      float alpha, 
                      float *A, magma_int_t lda, 
                      float *x, magma_int_t incx,
                      float beta,
                      float *z, magma_int_t incz)
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
             On entry, m specifies the number of rows of the matrix A.

    N      - (input) INTEGER
             On entry, n specifies the number of columns of the matrix A
 
    ALPHA  - REAL
             On entry, ALPHA specifies the scalar alpha.
             Unchanged on exit.

    A      - (input) SINGLE PRECISION array of dimension ( LDA, n ) on the GPU.
   
    LDA    - (input) INTEGER
             LDA specifies the leading dimension of A.

    X      - (input) SINGLE PRECISION array of dimension 
             n if trans == 'n'
             m if trans == 't'
     
    INCX   - (input) Specifies the increment for the elements of X.
             INCX must not be zero. Unchanged on exit.
  
    BETA   - REAL
             On entry, BETA specifies the scalar beta. When BETA is
             supplied as zero then Y need not be set on input.
             Unchanged on exit.

    Z      - (output) SINGLE PRECISION array of	dimension 
             m if trans == 'n'
             n if trans == 't' 

    INCZ  - (input) Specifies the increment for the elements of Z.
            INCZ must not be zero. Unchanged on exit.
    ===================================================================== */

    if (incx == 1 && incz == 1 && beta == 0.) {
       if (trans == 'n' || trans == 'N')
           magmablas_sgemvn_fermi(m,  n, alpha, A, lda, x, z);
       else if (trans == 't' || trans == 'T')
          magmablas_sgemvt_fermi(m,  n, alpha, A, lda, x, z);
       else
          printf("trans = %c in sgemv_fermi is not available\n", trans);	       
    }
    else
       cublasSgemv(trans, m, n, alpha, A, lda, x, incx, beta, z, incz);
}

#undef num_threads
#undef sgemv_bs
#undef threadSize 
