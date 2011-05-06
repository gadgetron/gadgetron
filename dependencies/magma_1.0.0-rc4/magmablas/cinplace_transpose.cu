/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
	Univ. of Colorado, Denver
       November 2010

       @generated c

*/
#include "common_magma.h"
#define PRECISION_c
#include "commonblas.h"

#if CSIZE_2SHARED == 32
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

__global__ void cinplace_T_even( cuFloatComplex *matrix, int lda, int half )
{	
        __shared__ cuFloatComplex a[CSIZE_2SHARED][CSIZE_2SHARED+1];
	__shared__ cuFloatComplex b[CSIZE_2SHARED][CSIZE_2SHARED+1];
	
	int inx = threadIdx.x;
	int iny = threadIdx.y;

	bool bottom = ( blockIdx.x > blockIdx.y );
	int ibx = bottom ? (blockIdx.x - 1) : (blockIdx.y + half);
	int iby = bottom ? blockIdx.y       : (blockIdx.x + half);

	ibx *= CSIZE_2SHARED;
	iby *= CSIZE_2SHARED;

	cuFloatComplex *A = matrix + ibx + inx + __mul24( iby + iny, lda );
	COPY_1D_TO_2D( A, lda, a, inx, iny);

	if( ibx == iby )
	{
		__syncthreads();
		COPY_2D_TO_1D( a, inx, iny, A, lda);
	}
	else
	{
		cuFloatComplex *B = matrix + iby + inx + __mul24( ibx + iny, lda );

		
		COPY_1D_TO_2D( B, lda, b, inx, iny);
		__syncthreads();

		COPY_2D_TO_1D( b, inx, iny, A, lda);
		COPY_2D_TO_1D( a, inx, iny, B, lda);
	}
} 

__global__ void cinplace_T_odd( cuFloatComplex *matrix, int lda, int half )
{	
        __shared__ cuFloatComplex a[CSIZE_2SHARED][CSIZE_2SHARED+1];
	__shared__ cuFloatComplex b[CSIZE_2SHARED][CSIZE_2SHARED+1];
	
	int inx = threadIdx.x;
	int iny = threadIdx.y;

	bool bottom = ( blockIdx.x >= blockIdx.y );
	int ibx = bottom ? blockIdx.x  : (blockIdx.y + half - 1);
	int iby = bottom ? blockIdx.y  : (blockIdx.x + half);

	ibx *= CSIZE_2SHARED;
	iby *= CSIZE_2SHARED;

	cuFloatComplex *A = matrix + ibx + inx + __mul24( iby + iny, lda );
	COPY_1D_TO_2D( A, lda, a, inx, iny);

	if( ibx == iby )
	{
		__syncthreads();
		COPY_2D_TO_1D( a, inx, iny, A, lda);
	}
	else
	{
		cuFloatComplex *B = matrix + iby + inx + __mul24( ibx + iny, lda );

		COPY_1D_TO_2D( B, lda, b, inx, iny);
		__syncthreads();

		COPY_2D_TO_1D( b, inx, iny, A, lda);
		COPY_2D_TO_1D( a, inx, iny, B, lda);
	}
} 

extern "C" void 
magmablas_cinplace_transpose( cuFloatComplex *A, int lda, int n )
{
	dim3 threads( CSIZE_2SHARED, 16 );
	int in = n / CSIZE_2SHARED;
	if( in&1 )
	{
		dim3 grid( in, in/2+1 );
		cinplace_T_odd<<< grid, threads >>>( A, lda, in/2+1 );
	}
	else
	{
		dim3 grid( in+1, in/2 );
		cinplace_T_even<<< grid, threads >>>( A, lda, in/2 );
	}
}
