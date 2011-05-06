/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"
#include "commonblas_d.h"

static __device__ void daxpy(double a,double *b, double *c) {
	c[0] += a * b[0];
	c[1] += a * b[1];
	c[2] += a * b[2];
	c[3] += a * b[3];
	c[4] += a * b[4];
	c[5] += a * b[5];
	c[6] += a * b[6];
	c[7] += a * b[7];
	c[8] += a * b[8];
	c[9] += a * b[9];
	c[10] += a * b[10];
	c[11] += a * b[11];
	c[12] += a * b[12];
	c[13] += a * b[13];
	c[14] += a * b[14];
	c[15] += a * b[15];
}

extern "C" __global__ void 
dgemm_kernel_N_N_64_16_16_16_4_special(double *C, const double *A, const double *B,
                                       int m, int n, int k, 
                                       int lda, int ldb, int ldc, 
                                       double alpha, double beta) 
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose:
    ========
    This routine computes
       C = alpha* A*B  + beta * C

    B is put into shared memory
    Parameters Used:
        blk_M=64 blk_N=16 blk_K=16 nthd_x=16 nthd_y=4

    This kernel is for matrices devisible by the corresponding
    blocking sizes.
    ===============================================================  */

	const int tx = threadIdx.x;
	const int ty = threadIdx.y;

	const int ibx = blockIdx.x * 64;
	const int iby = blockIdx.y *16;

	const int idt = ty * 16 + tx;

	B+=tx+__mul24(iby+ty,ldb);
	A += ibx + idt;
	C += ibx +idt +__mul24( iby,ldc);

	const double *Bend = B + k;


	double Cb[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	m = 2*lda ; 
	n = 3*lda ;

	do {
		//double Ab[4] = {A[0], A[lda], A[2*lda], A[3*lda]};
		double Ab[4] = {A[0], A[lda], A[m], A[n]};
		__shared__ double Bb[16][17];
		Bb[tx][ty+0] = B[0];
		Bb[tx][ty+4] = B[4*ldb];
		Bb[tx][ty+8] = B[8*ldb];
		Bb[tx][ty+12] = B[12*ldb];

		__syncthreads();

		A += 4 * lda;
		daxpy(Ab[0], &Bb[0][0], Cb); Ab[0] = A[0];
		daxpy(Ab[1], &Bb[1][0], Cb); Ab[1] = A[lda];
		daxpy(Ab[2], &Bb[2][0], Cb); Ab[2] = A[m];
		daxpy(Ab[3], &Bb[3][0], Cb); Ab[3] = A[n];

		A += 4 * lda;
		daxpy(Ab[0], &Bb[4][0], Cb); Ab[0] = A[0];
		daxpy(Ab[1], &Bb[5][0], Cb); Ab[1] = A[lda];
		daxpy(Ab[2], &Bb[6][0], Cb); Ab[2] = A[m];
		daxpy(Ab[3], &Bb[7][0], Cb); Ab[3] = A[n];

		A += 4 * lda;
		daxpy(Ab[0], &Bb[8][0], Cb); Ab[0] = A[0];
		daxpy(Ab[1], &Bb[9][0], Cb); Ab[1] = A[lda];
		daxpy(Ab[2], &Bb[10][0], Cb); Ab[2] = A[m];
		daxpy(Ab[3], &Bb[11][0], Cb); Ab[3] = A[n];

		A += 4 * lda;
		daxpy(Ab[0], &Bb[12][0], Cb);
		daxpy(Ab[1], &Bb[13][0], Cb);
		daxpy(Ab[2], &Bb[14][0], Cb);
		daxpy(Ab[3], &Bb[15][0], Cb);

		B += 16;

		__syncthreads();
	} while (B < Bend);

	#pragma unroll 16
	for (int i = 0; i < 16; i++, C += ldc) {
		C[0] =alpha*Cb[i] + beta * C[0];
	}
}

extern "C" void
magmablas_dgemm_kernel_N_N_64_16_16_16_4_special(double *C, 
                                                 const double *A, 
                                                 const double *B, 
                                                 int m, int n, int k,
                                                 int lda, int ldb, int ldc,
                                                 double alpha, double beta)
{
        dim3 threads( 16, 4 );
        dim3 grid(m/64,n/16);
        dgemm_kernel_N_N_64_16_16_16_4_special<<< grid, threads >>>(C, A, B,
                                                                    m, n, k,
                                                                    lda, ldb, ldc, 
                                                                    alpha, beta);
}

