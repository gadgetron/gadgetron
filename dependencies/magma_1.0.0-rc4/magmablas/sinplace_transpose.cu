/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
	Univ. of Colorado, Denver
       November 2010

       @generated s

*/
#include "common_magma.h"
#define PRECISION_s
#include "commonblas.h"

#if SSIZE_2SHARED == 32
#define COPY_1D_TO_2D( a, lda, b, inx, iny )	\
    b[iny][inx]    = a[0];			\
    b[iny+16][inx] = a[16*lda];

#define COPY_2D_TO_1D( b, inx, iny, a, lda )	\
    a[0]      = b[inx][iny];			\
    a[16*lda] = b[inx][iny+16];

#else

#define COPY_1D_TO_2D( a, lda, b, inx, iny )	\
    b[iny][inx] = a[0];

#define COPY_2D_TO_1D( b, inx, iny, a, lda )	\
    a[0]        = b[inx][iny];

#endif

__global__ void sinplace_T_even( float *matrix, int lda, int half )
{	
        __shared__ float a[SSIZE_2SHARED][SSIZE_2SHARED+1];
	__shared__ float b[SSIZE_2SHARED][SSIZE_2SHARED+1];
	
	int inx = threadIdx.x;
	int iny = threadIdx.y;

	bool bottom = ( blockIdx.x > blockIdx.y );
	int ibx = bottom ? (blockIdx.x - 1) : (blockIdx.y + half);
	int iby = bottom ? blockIdx.y       : (blockIdx.x + half);

	ibx *= SSIZE_2SHARED;
	iby *= SSIZE_2SHARED;

	float *A = matrix + ibx + inx + __mul24( iby + iny, lda );
	COPY_1D_TO_2D( A, lda, a, inx, iny);

	if( ibx == iby )
	{
		__syncthreads();
		COPY_2D_TO_1D( a, inx, iny, A, lda);
	}
	else
	{
		float *B = matrix + iby + inx + __mul24( ibx + iny, lda );

		
		COPY_1D_TO_2D( B, lda, b, inx, iny);
		__syncthreads();

		COPY_2D_TO_1D( b, inx, iny, A, lda);
		COPY_2D_TO_1D( a, inx, iny, B, lda);
	}
} 

__global__ void sinplace_T_odd( float *matrix, int lda, int half )
{	
        __shared__ float a[SSIZE_2SHARED][SSIZE_2SHARED+1];
	__shared__ float b[SSIZE_2SHARED][SSIZE_2SHARED+1];
	
	int inx = threadIdx.x;
	int iny = threadIdx.y;

	bool bottom = ( blockIdx.x >= blockIdx.y );
	int ibx = bottom ? blockIdx.x  : (blockIdx.y + half - 1);
	int iby = bottom ? blockIdx.y  : (blockIdx.x + half);

	ibx *= SSIZE_2SHARED;
	iby *= SSIZE_2SHARED;

	float *A = matrix + ibx + inx + __mul24( iby + iny, lda );
	COPY_1D_TO_2D( A, lda, a, inx, iny);

	if( ibx == iby )
	{
		__syncthreads();
		COPY_2D_TO_1D( a, inx, iny, A, lda);
	}
	else
	{
		float *B = matrix + iby + inx + __mul24( ibx + iny, lda );

		COPY_1D_TO_2D( B, lda, b, inx, iny);
		__syncthreads();

		COPY_2D_TO_1D( b, inx, iny, A, lda);
		COPY_2D_TO_1D( a, inx, iny, B, lda);
	}
} 

extern "C" void 
magmablas_sinplace_transpose( float *A, int lda, int n )
{
	dim3 threads( SSIZE_2SHARED, 16 );
	int in = n / SSIZE_2SHARED;
	if( in&1 )
	{
		dim3 grid( in, in/2+1 );
		sinplace_T_odd<<< grid, threads >>>( A, lda, in/2+1 );
	}
	else
	{
		dim3 grid( in+1, in/2 );
		sinplace_T_even<<< grid, threads >>>( A, lda, in/2 );
	}
}
