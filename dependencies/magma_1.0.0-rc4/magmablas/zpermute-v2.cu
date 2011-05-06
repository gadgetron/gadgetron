/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        cuDoubleComplex *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} zlaswp_params_t;

typedef struct {
        cuDoubleComplex *A;
        int n, lda, j0, npivots;
        short ipiv[BLOCK_SIZE];
} zlaswp_params_t2;

/*********************************************************
 *
 * LAPACK Swap: permute a set of lines following ipiv
 *
 ********************************************************/
typedef struct {
    cuDoubleComplex *A;
    int n, ldx, ldy, j0, npivots;
    short ipiv[BLOCK_SIZE];
} zlaswpx_params_t;

__global__ void myzlaswpx( zlaswpx_params_t params )
{
    unsigned int y = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
    unsigned int offset1 = __mul24( y, params.ldy);
    if( y < params.n )
    {
        int ldx = params.ldx;
        cuDoubleComplex *A = params.A + offset1 + ldx * params.j0;
        cuDoubleComplex *Ai = A;
        
        for( int i = 0; i < params.npivots; i++ )
        {
            int j = params.ipiv[i];
            cuDoubleComplex *p2 = A + j*ldx;
            cuDoubleComplex temp = *Ai;
            *Ai = *p2;
            *p2 = temp;
            Ai += ldx;
        }
    }
}

extern "C" void zlaswpx( zlaswpx_params_t &params )
{
 	int blocksize = 64;
	dim3 blocks = (params.n+blocksize-1) / blocksize;
	myzlaswpx<<< blocks, blocksize >>>( params );
}

/*
 * Old version
 */
__global__ void myzlaswp2( zlaswp_params_t2 params )
{
        unsigned int tid = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
        if( tid < params.n )
	{
                int lda = params.lda;
		cuDoubleComplex *A = params.A + tid + lda * params.j0;

		for( int i = 0; i < params.npivots; i++ )
		{
                 	int j = params.ipiv[i];
			cuDoubleComplex *p1 = A + i*lda;
			cuDoubleComplex *p2 = A + j*lda;
			cuDoubleComplex temp = *p1;
			*p1 = *p2;
			*p2 = temp;
		}
	}
}

extern "C" void zlaswp2( zlaswp_params_t &params );

extern "C" void zlaswp3( zlaswp_params_t2 &params )
{
 	int blocksize = 64;
	dim3 blocks = (params.n+blocksize-1) / blocksize;
	myzlaswp2<<< blocks, blocksize >>>( params );
}


extern "C" void 
magmablas_zpermute_long2( cuDoubleComplex *dAT, int lda, int *ipiv, int nb, int ind )
{
        int k;

        for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
        {
                //zlaswp_params_t params = { dAT, lda, lda, ind + k };
                zlaswp_params_t2 params = { dAT, lda, lda, ind + k, BLOCK_SIZE };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1;
                        ipiv[ind + k + j] += ind;
                }
                //zlaswp2( params );
	        zlaswp3( params );
        }

	int num_pivots = nb - k;

        zlaswp_params_t2 params = { dAT, lda, lda, ind + k, num_pivots};
        for( int j = 0; j < num_pivots; j++ )
        {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1;
            ipiv[ind + k + j] += ind;
        }
        zlaswp3( params );
}

extern "C" void 
magmablas_zlaswp( int n, cuDoubleComplex *dAT, int lda, 
                  int i1, int i2, int *ipiv, int inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
    {
      int sb = min(BLOCK_SIZE, i2-k);
      //zlaswp_params_t params = { dAT, lda, lda, ind + k };
      zlaswp_params_t2 params = { dAT+k*lda, n, lda, 0, sb };
      for( int j = 0; j < sb; j++ )
        {
          params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
      zlaswp3( params );
    }
}

extern "C" void 
magmablas_zlaswpx( int n, cuDoubleComplex *dAT, int ldx, int ldy, 
                   int i1, int i2, int *ipiv, int inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
    {
      int sb = min(BLOCK_SIZE, i2-k);
      //zlaswp_params_t params = { dAT, lda, lda, ind + k };
      zlaswpx_params_t params = { dAT+k*ldx, n, ldx, ldy, 0, sb };
      for( int j = 0; j < sb; j++ )
        {
          params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
      zlaswpx( params );
    }
}

#undef BLOCK_SIZE
