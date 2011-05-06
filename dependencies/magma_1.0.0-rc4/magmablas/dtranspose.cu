/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d

*/
#include "common_magma.h"
#define PRECISION_d
#include "commonblas.h"

__global__ void dtranspose_32( double *B, int ldb, double *A, int lda )
{	
	__shared__ double a[32][DSIZE_1SHARED+1];
	
	int inx = threadIdx.x;
	int iny = threadIdx.y;
	int ibx = blockIdx.x*32;
	int iby = blockIdx.y*32;
	
	A += ibx + inx + __mul24( iby + iny, lda );
	B += iby + inx + __mul24( ibx + iny, ldb );
	
	a[iny+0][inx] = A[0*lda];
	a[iny+8][inx] = A[8*lda];
	a[iny+16][inx] = A[16*lda];
	a[iny+24][inx] = A[24*lda];
	
	__syncthreads();
	
#if defined(PRECISION_s) || defined(PRECISION_d) || defined(PRECISION_c)
	B[0*ldb] = a[inx][iny+0];
	B[8*ldb] = a[inx][iny+8];
	B[16*ldb] = a[inx][iny+16];
	B[24*ldb] = a[inx][iny+24];
#else /* defined(PRECISION_z) */
	B[0*ldb]    = a[inx][iny+0];
	B[8*ldb]    = a[inx][iny+8];
	B[0*ldb+16] = a[inx+16][iny+0];
	B[8*ldb+16] = a[inx+16][iny+8];

	__syncthreads();
	A += DSIZE_1SHARED;
	B += __mul24( 16, ldb);

        a[iny+0][inx] = A[0*lda];
        a[iny+8][inx] = A[8*lda];
        a[iny+16][inx] = A[16*lda];
        a[iny+24][inx] = A[24*lda];

        __syncthreads();

	B[0*ldb] = a[inx][iny+0];
	B[8*ldb] = a[inx][iny+8];
	B[0*ldb+16] = a[inx+16][iny+0];
	B[8*ldb+16] = a[inx+16][iny+8];
#endif
} 

//
//	m, n - dimensions in the source matrix
//             This version works when m and n are divisible by 32.
//
extern "C" void 
magmablas_dtranspose(double *odata, int ldo, 
                     double *idata, int ldi, 
                     int m, int n )
{
	//assert( (m%32) == 0 && (n%32) == 0, "misaligned transpose" );
	dim3 threads( DSIZE_1SHARED, 8, 1 );
	dim3 grid( m/32, n/32, 1 );
	dtranspose_32<<< grid, threads >>>( odata, ldo, idata, ldi );
}
