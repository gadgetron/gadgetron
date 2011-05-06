/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"

#define num_threads 128
#define zgemv_bs     32
#define threadSize  128

#define magmablas_zgemv_fermi magmablas_zgemv

__global__ void 
zgemvn_kernel1_fermi(int n, int m, int n1, cuDoubleComplex alpha, cuDoubleComplex* A, int lda, cuDoubleComplex *x, cuDoubleComplex *y)
{
  int ind = blockIdx.x*num_threads + threadIdx.x;

  A += ind;

  cuDoubleComplex res;
  MAGMA_Z_SET2REAL(res, 0.0f);

  for(int i=0; i<n1; i += zgemv_bs ){

    #pragma unroll
    for(int j=0; j < zgemv_bs ; j++){
       res += A[0] * x[j];
       A   += lda;
    }
	x += zgemv_bs;
  }

  if (m>n1){

     for(int j=0; j<(m-n1); j++){
         res += A[0] * x[j];
         A   += lda;
     }
  }

  if (ind<n)
     y[ind] = alpha * res;

}

__global__ void 
zgemvn_kernel2_fermi(int n, int m, int n1, cuDoubleComplex alpha,  cuDoubleComplex* A, int lda, cuDoubleComplex *x, cuDoubleComplex *y)
{
  int ind = blockIdx.x*num_threads + threadIdx.x;

  A += ind;
  x += threadIdx.x;

  cuDoubleComplex res;
  MAGMA_Z_SET2REAL(res, 0.0f);

  __shared__ cuDoubleComplex buff[num_threads];
  for(int i=0; i<n1; i += num_threads ){
    __syncthreads();
    buff[threadIdx.x]  = x[i];

    __syncthreads();
    #pragma unroll
    for(int j=0; j < num_threads ; j++){
       res+=A[0]*buff[j];
       A+=lda;
    }
  }
  __syncthreads();

  if (m>n1){
     buff[threadIdx.x]  = x[n1];

     __syncthreads();
     for(int j=0; j<(m-n1); j++){
         res += A[0]*buff[j];
         A+=lda;
     }
  }

  if (ind<n)
     y[ind] = alpha * res;
}

extern "C" void
magmablas_zgemvn_fermi(int n, int m, cuDoubleComplex alpha, cuDoubleComplex *A, int lda, cuDoubleComplex *x, cuDoubleComplex *y)
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

    int blocks;
    if (n % num_threads==0)
        blocks = n/num_threads;
    else
        blocks = n/num_threads + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(num_threads, 1, 1);
/*    if(n<=8500) 
		zgemvn_kernel1_fermi<<<grid, threads>>>(n, m, (m / zgemv_bs)*zgemv_bs, 
			                           alpha, A, lda, x, y);
    else */
		zgemvn_kernel2_fermi<<<grid, threads>>>(n, m, (m / num_threads)*num_threads, 
			                           alpha, A, lda, x, y);
}



__global__ void 
zgemvt_kernel_fermi(int m, int n, cuDoubleComplex alpha, int n1, cuDoubleComplex* A, int lda,
              cuDoubleComplex *x, cuDoubleComplex *y)
{
	unsigned int tx = threadIdx.x;

	__shared__ cuDoubleComplex sdata[threadSize];
	

	cuDoubleComplex res;
    MAGMA_Z_SET2REAL(res, 0.0);
	cuDoubleComplex zero;
    MAGMA_Z_SET2REAL(zero, 0.0);
     
	for(int i=0; i<n1; i+= threadSize)
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
			res  += zero;
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
		sdata[tx] += sdata[tx + 32];
	}

    if(tx == 0)
	{
		for(int i=1;i<32;i++)
		{
			sdata[tx] += sdata[tx + i];
		}
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




extern "C" void
magmablas_zgemvt_fermi(int m, int n, cuDoubleComplex alpha, cuDoubleComplex *A, int lda, 
                 cuDoubleComplex *x, cuDoubleComplex *y)
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

    
    dim3 grid    ( 1,  n,  1);
    dim3 threads ( threadSize,   1,  1);

    zgemvt_kernel_fermi<<<grid, threads>>>( m, n, alpha, ( m / threadSize) * threadSize,
                                       A, lda, x, y);

}



__global__ void 
zgemvc_kernel_fermi(int m, int n, cuDoubleComplex alpha, int n1, cuDoubleComplex* A, int lda,
              cuDoubleComplex *x, cuDoubleComplex *y)
{
	unsigned int tx = threadIdx.x;

	__shared__ cuDoubleComplex sdata[threadSize];
	

	cuDoubleComplex res;
    MAGMA_Z_SET2REAL(res, 0.0);
	cuDoubleComplex zero;
    MAGMA_Z_SET2REAL(zero, 0.0);
     
	for(int i=0; i<n1; i+= threadSize)
	{
		res += cuConj(A[tx + i + lda * blockIdx.y]) * x[tx + i];
	}

	
	if(m > n1)
	{
		if( tx + n1 <  m )
		{
			res  += cuConj(A[tx + n1 + lda *blockIdx.y]) * x[tx + n1];
		}
		else 
		{
			res  += zero;
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
		sdata[tx] += sdata[tx + 32];
	}

    if(tx == 0)
	{
		for(int i=1;i<32;i++)
		{
			sdata[tx] += sdata[tx + i];
		}
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




extern "C" void
magmablas_zgemvc_fermi(int m, int n, cuDoubleComplex alpha, cuDoubleComplex *A, int lda, 
                 cuDoubleComplex *x, cuDoubleComplex *y)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    This routine computes y = alpha *conjg(A^t) *  x on the GPU.

    M      - (input) INTEGER.
             On entry, M specifies the number of rows of the matrix A.

    N      - (input) INTEGER.
             On entry, N specifies the number of columns of the matrix A

    A      - (input) SINGLE PRECISION array of dimension ( LDA, n ) on the GPU.

    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) SINGLE PRECISION array of dimension m.

    Y      - (output) SINGLE PRECISION array of dimension n.
             On exit Y = alpha conjg(A^t) X.

    ===================================================================== */

    
    dim3 grid    ( 1,  n,  1);
    dim3 threads ( threadSize,   1,  1);

    zgemvc_kernel_fermi<<<grid, threads>>>( m, n, alpha, ( m / threadSize) * threadSize,
                                       A, lda, x, y);

}





extern "C" void
magmablas_zgemv_fermi(char flag, int m, int n, cuDoubleComplex alpha, cuDoubleComplex *A, int lda, cuDoubleComplex *x, int incx, cuDoubleComplex beta, cuDoubleComplex *y, int incy) 
{

    if(beta.x==0.0 && beta.y==0.0 && incx ==1 && incy==1)
	{
		if (flag == 'N' || flag == 'n')
		{
			if(m<7000)
			{
				cublasZgemv(flag, m, n, alpha, A, lda, x, incx, beta, y, incy);
			}
			else
			{
				magmablas_zgemvn_fermi(m,  n, alpha, A, lda, x, y);
			}
		}
		else if(flag == 'T' || flag == 't')
		{
			magmablas_zgemvt_fermi(m,  n, alpha, A, lda, x, y);
		}
		else if(flag == 'C' || flag == 'c') 
		{
			magmablas_zgemvc_fermi(m,  n, alpha, A, lda, x, y);
		}
		else 
		{
			cublasZgemv(flag, m, n, alpha, A, lda, x, incx, beta, y, incy);
		}
	}
	else
	{
		cublasZgemv(flag, m, n, alpha, A, lda, x, incx, beta, y, incy);
	}
}

#undef num_threads
#undef zgemv_bs
#undef threadSize 
