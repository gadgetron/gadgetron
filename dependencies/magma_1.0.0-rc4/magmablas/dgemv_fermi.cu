
/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"

#define magmablas_dgemv_fermi magmablas_dgemv

#define num_threads 64
#define dgemv_bs 64
#define threadSize 128



__global__ void 
dgemvn_kernel_fermi(magma_int_t n, magma_int_t m, magma_int_t n1, double alpha,  double* A, magma_int_t lda, double *x, double *y)
{
  magma_int_t ind = blockIdx.x*num_threads + threadIdx.x;

  A += ind;
  x += threadIdx.x;


  double res = 0.0;

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
     y[ind] = alpha * res;
}


extern "C" void
magmablas_dgemvn_fermi(magma_int_t n, magma_int_t m, double alpha, double *A, magma_int_t lda, double *x, double *y)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    This routine computes y = alpha A x on the GPU.

    N      - (input) INTEGER.
             On entry, N specifies the number of rows of the matrix A.

    M      - (input) INTEGER.
             On entry, M specifies the number of columns of the matrix A

    A      - (input) DOUBLE PRECISION array of dimension ( LDA, m ) on the GPU.
   
    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) DOUBLE PRECISION array of dimension m.
     
    Y      - (output) DOUBLE PRECISION array of	dimension m. 
             On exit Y = alpha A X.

    ===================================================================== */

    magma_int_t blocks;
    if (n % num_threads==0)
        blocks = n/num_threads;
    else
        blocks = n/num_threads + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(num_threads, 1, 1);
 
    dgemvn_kernel_fermi<<<grid, threads>>>(n, m, (m / dgemv_bs)*dgemv_bs, 
                                    alpha, A, lda, x, y);
}


__global__ void
dgemvt_kernel_fermi(magma_int_t m, magma_int_t n, double alpha, magma_int_t n1, double* A, magma_int_t lda,
              double *x, double *y)
{
	magma_int_t tx = threadIdx.x;

	__shared__ double sdata[threadSize];


	double res;
	res = 0.0;

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
			res  = res;
		}
    }

	sdata[tx] = res;
	 __syncthreads();


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
magmablas_dgemvt_fermi(magma_int_t m, magma_int_t n, double alpha, double *A, magma_int_t lda,
                 double *x, double *y)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    This routine computes y = alpha A^t x on the GPU.

    M      - (input) INTEGER.
             On entry, m specifies the number of rows of the matrix A.

    N      - (input) INTEGER.
             On entry, n specifies the number of columns of the matrix A

    A      - (input) DOUBLE PRECISION array of dimension ( LDA, n ) on the GPU.

    LDA    - (input) INTEGER.
             LDA specifies the leading dimension of A.

    X      - (input) DOUBLE PRECISION array of dimension m.

    Y      - (output) DOUBLE PRECISION array of dimension n.
             On exit y = alpha A^t X.

    ===================================================================== */

	dim3 grid    ( 1,  n,  1);
	dim3 threads ( threadSize,   1,  1);

	dgemvt_kernel_fermi<<<grid, threads>>>( m, n, alpha, ( m / threadSize) * threadSize,
				                                       A, lda, x, y);
}



extern "C" 
void magmablas_dgemv_fermi(char trans,
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
             On entry, m specifies the number of rows of the matrix A.

    N      - (input) INTEGER
             On entry, n specifies the number of columns of the matrix A
 
    ALPHA  - DOUBLE REAL
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
  
    BETA   - DOUBLE REAL
             On entry, BETA specifies the scalar beta. When BETA is
             supplied as zero then Y need not be set on input.
             Unchanged on exit.

    Z      - (output) DOUBLE PRECISION array of	dimension 
             m if trans == 'n'
             n if trans == 't' 

    INCZ  - (input) Specifies the increment for the elements of Z.
            INCZ must not be zero. Unchanged on exit.
    ===================================================================== */

    if (incx == 1 && incz == 1 && beta == 0.) {
       if (trans == 'n' || trans == 'N')
	   {
	       if ( m >= 7000 && m <= 8000 )
                cublasDgemv(trans, m, n, alpha, A, lda, x, incx, beta, z, incz);
		   else 
				magmablas_dgemvn_fermi(m,  n, alpha, A, lda, x, z);
	   }
       else if (trans == 't' || trans == 'T')
          magmablas_dgemvt_fermi(m,  n, alpha, A, lda, x, z);
       else
          printf("trans = %c in sgemv_fermi is not available\n", trans);	       
    }
    else
       cublasDgemv(trans, m, n, alpha, A, lda, x, incx, beta, z, incz);
}

#undef num_threads
#undef dgemv_bs
#undef threadSize
