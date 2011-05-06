/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated c

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

/*********************************************************/
/*
*  Blocked version: swap several pair of line
 */
typedef struct {
    cuFloatComplex *A1;
    cuFloatComplex *A2;
    int n, lda1, lda2, npivots;
    short ipiv[BLOCK_SIZE];
} magmagpu_cswapblk_params_t;

__global__ void magmagpu_cswapblkrm( magmagpu_cswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    if( y < params.n )
    {
        cuFloatComplex *A1 = params.A1 + y - params.lda1;
        cuFloatComplex *A2 = params.A2 + y;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A1 += params.lda1;
            if ( params.ipiv[i] == -1 )
                continue;
            cuFloatComplex tmp1  = *A1;
            cuFloatComplex *tmp2 = A2 + params.ipiv[i]*params.lda2;
            *A1   = *tmp2;
            *tmp2 = tmp1;
        }
    }
}

__global__ void magmagpu_cswapblkcm( magmagpu_cswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned int offset1 = __mul24( y, params.lda1);
    unsigned int offset2 = __mul24( y, params.lda2);
    if( y < params.n )
    {
        cuFloatComplex *A1 = params.A1 + offset1 - 1;
        cuFloatComplex *A2 = params.A2 + offset2;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A1++;
            if ( params.ipiv[i] == -1 )
                continue;
            cuFloatComplex tmp1  = *A1;
            cuFloatComplex *tmp2 = A2 + params.ipiv[i];
            *A1   = *tmp2;
            *tmp2 = tmp1;
        }
    }
    __syncthreads();
}

extern "C" void 
magmablas_cswapblk( char storev, int n, 
                    cuFloatComplex *dA1T, int lda1,
                    cuFloatComplex *dA2T, int lda2,
                    int i1, int i2, int *ipiv, int inci, int offset )
{
    int  blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    int  k, im;

    if ( (storev == 'C') || (storev == 'c') ) {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            int sb = min(BLOCK_SIZE, i2-k);
            magmagpu_cswapblk_params_t params = { dA1T+k, dA2T, n, lda1, lda2, sb };
            for( int j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im)
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magmagpu_cswapblkcm<<< blocks, blocksize >>>( params );
        }
    }else {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            int sb = min(BLOCK_SIZE, i2-k);
            magmagpu_cswapblk_params_t params = { dA1T+k*lda1, dA2T, n, lda1, lda2, sb };
            for( int j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im)
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magmagpu_cswapblkrm<<< blocks, blocksize >>>( params );
        }
    }
}

