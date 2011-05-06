/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       February 2011
*/
#include "common_magma.h"

#define magmablas_dtrsm_tesla magmablas_dtrsm

#define qmod(a,b) ((a)-(__mul24((b),(a)/(b))))
#define cublasDgemm magmablas_dgemm

#define b_copy();	dim3 dimBlock((M>=MAX_THREAD_PER_BLOCK)?MAX_THREAD_PER_BLOCK:(WARP_SIZE*((M/WARP_SIZE)+(M%WARP_SIZE!=0))), 1);\
					dim3 dimGrid(M/dimBlock.x+(M%dimBlock.x!=0), N);\
					b_copy_kernel<<<dimGrid, dimBlock>>>(M, N, b, ldb, d_x, M);\
					cudaThreadSynchronize();

#define MAX_THREAD_PER_BLOCK 512
#define WARP_SIZE 32

#define BLOCK_SIZE 16 // inner blocking size, <=32
#define NB 128// outer blocking size, >BLOCK_SIZE 

__global__ void
diag_dtrtri_kernel_upper (char diag, double *A, double *d_dinvA, int lda)
{
	int i,j;
	double Ystx=0;
	double *y=NULL, *Aoff=NULL;
	int switcher=0;

	// Thread index
	int tx = threadIdx.x;

	// Block index
	int bx = blockIdx.x;
		
	Aoff = A+bx*lda*BLOCK_SIZE+bx*BLOCK_SIZE;
	int NumBLperNB = NB/BLOCK_SIZE;
	d_dinvA += bx/NumBLperNB*NB*NB+(bx%NumBLperNB)*(NB*BLOCK_SIZE+BLOCK_SIZE);

	__shared__ double Bs[BLOCK_SIZE*BLOCK_SIZE];
	__shared__ double workspace[BLOCK_SIZE]; // workspace used to store the current working column

	// load A
	#pragma unroll
	for (i=0; i<BLOCK_SIZE; i++)
		Bs[i*BLOCK_SIZE+tx] = ((double)(tx<=i))*(*(Aoff+i*lda+tx));	// read in the whole square block of my A and zero out the non data triangular

	// Synchronize to make sure the matrices are loaded
	__syncthreads();

	switcher = (diag=='u' || diag=='U');
	int diagsw = (Bs[tx*BLOCK_SIZE+tx]==0);
	Bs[tx*BLOCK_SIZE+tx] = switcher+!switcher*(1/(diagsw+(!diagsw)*Bs[tx*BLOCK_SIZE+tx]));	// solve the diagonals

	/* the upper case */
	for (i=0; i<BLOCK_SIZE; i++)
	{
		Ystx = 0;
		switcher = (double)(tx<i);

		//dtrmv
		workspace[tx] = *(Bs+i*BLOCK_SIZE+tx);
		y = Bs+i*BLOCK_SIZE;

		#pragma unroll
		//for (j=tx; j<i; j++)
		for (j=0; j<i; j++)
			Ystx += switcher*(*(Bs+j*BLOCK_SIZE+tx)*workspace[j]);

		//sscal
		switcher = (tx != i); // if (tx !=i ) y[tx]=switcher*Ystx*(-Bs[i*BLOCK_SIZE+i]);
		y[tx] = switcher*Ystx*(-Bs[i*BLOCK_SIZE+i])+!switcher*y[tx];

		__syncthreads();
	}

	// write back A
	#pragma unroll
	for (i=0; i<BLOCK_SIZE; i++)
		*(d_dinvA+i*NB+tx) = Bs[i*BLOCK_SIZE+tx];

}
__global__ void
diag_dtrtri_kernel_lower (char diag, double *A, double *d_dinvA, int lda)
{
	int i,j;
	double Ystx=0;
	double *Bw=NULL, *x=NULL, *y=NULL, *Aoff=NULL;
	int switcher=0;

	// Thread index
	int tx = threadIdx.x;
	int txw;

	// Block index
	int bx = blockIdx.x;
		
	Aoff = A+bx*lda*BLOCK_SIZE+bx*BLOCK_SIZE;
	int NumBLperNB = NB/BLOCK_SIZE;
	d_dinvA += bx/NumBLperNB*NB*NB+(bx%NumBLperNB)*(NB*BLOCK_SIZE+BLOCK_SIZE);

	__shared__ double Bs[BLOCK_SIZE*BLOCK_SIZE];
	__shared__ double workspace[BLOCK_SIZE]; // workspace used to store the current working column

	// load A
	#pragma unroll
	for (i=0; i<BLOCK_SIZE; i++)
		Bs[i*BLOCK_SIZE+tx] = ((double)(tx>=i))*(*(Aoff+i*lda+tx));	// read in the whole square block of my A and zero out the non data triangular
												// not the upper or lower diagonal
	// Synchronize to make sure the matrices are loaded
	__syncthreads();

	switcher = (diag=='u' || diag=='U');
	int diagsw = (Bs[tx*BLOCK_SIZE+tx]==0);
	Bs[tx*BLOCK_SIZE+tx] = switcher+!switcher*(1/(diagsw+(!diagsw)*Bs[tx*BLOCK_SIZE+tx]));	// solve the diagonals

	/*
	 * the lower case
	 */

	switcher = !(tx<BLOCK_SIZE-1);
	Bs[(BLOCK_SIZE-1)*BLOCK_SIZE+tx] = (double)switcher*Bs[(BLOCK_SIZE-1)*BLOCK_SIZE+tx];	//zero out the last column, except the diagonal element

	for (i=BLOCK_SIZE-2; i>=0; i--)
	{
		Ystx = 0;
		switcher = (tx>i);

		//dtrmv
		Bw = Bs+(i+1)*BLOCK_SIZE+i+1;
		workspace[tx] = *(Bs+i*BLOCK_SIZE+tx);
		x = workspace+i+1;
		y = Bs+i*BLOCK_SIZE;

		txw = (tx-i-1);

		#pragma unroll
		for (j=0; j<BLOCK_SIZE-i-1; j++)
			Ystx += (double)switcher*(*(Bw+j*BLOCK_SIZE+txw)*x[j]);

		//sscal
		switcher = (tx != i); 
		y[tx] = (double)switcher*Ystx*(-Bs[i*BLOCK_SIZE+i])+(double)(!switcher)*y[tx];

		__syncthreads();
	}

	// write back A
	#pragma unroll
	for (i=0; i<BLOCK_SIZE; i++)
		*(d_dinvA+i*NB+tx) = Bs[i*BLOCK_SIZE+tx];
}

__device__ void ssaxpy( double a, double *b, double *c )
{
	c[0] += a*b[0];
	c[1] += a*b[1];
	c[2] += a*b[2];
	c[3] += a*b[3];
	c[4] += a*b[4];
	c[5] += a*b[5];
	c[6] += a*b[6];
	c[7] += a*b[7];
	c[8] += a*b[8];
	c[9] += a*b[9];
	c[10] += a*b[10];
	c[11] += a*b[11];
	c[12] += a*b[12];
	c[13] += a*b[13];
	c[14] += a*b[14];
	c[15] += a*b[15];
}

__device__ void dgemm_kernel_16 (double *A, int lda, double *B, int ldb, double * C, int ldc, double alpha, int blk, int inx, int iny, double *c)
{

	const double *Blast = B + blk;
	__shared__ double bs[16][17];


	do
	{
		double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

		bs[inx][iny]    = B[0*ldb];
		bs[inx][iny+4]  = B[4*ldb];
		bs[inx][iny+8]  = B[8*ldb];
		bs[inx][iny+12] = B[12*ldb];
		bs[inx+4][iny]    = B[4+0*ldb];
		bs[inx+4][iny+4]  = B[4+4*ldb];
		bs[inx+4][iny+8]  = B[4+8*ldb];
		bs[inx+4][iny+12] = B[4+12*ldb];
		bs[inx+8][iny]    = B[8+0*ldb];
		bs[inx+8][iny+4]  = B[8+4*ldb];
		bs[inx+8][iny+8]  = B[8+8*ldb];
		bs[inx+8][iny+12] = B[8+12*ldb];
		bs[inx+12][iny]    = B[12+0*ldb];
		bs[inx+12][iny+4]  = B[12+4*ldb];
		bs[inx+12][iny+8]  = B[12+8*ldb];
		bs[inx+12][iny+12] = B[12+12*ldb];
		__syncthreads();

		A += 4*lda;
		ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
		ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
		ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
		ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

		A += 4*lda;
		ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
		ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
		ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
		ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

		A += 4*lda;
		ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
		ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
		ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
		ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

		A += 4*lda;
		ssaxpy( a[0], &bs[12][0], c );
		ssaxpy( a[1], &bs[13][0], c );
		ssaxpy( a[2], &bs[14][0], c );
		ssaxpy( a[3], &bs[15][0], c );

		B += 16;
		__syncthreads();
	} while( B < Blast );

	for( int i = 0; i < 16; i++, C += ldc )
		C[0] = alpha*c[i];
}

/*
 * B21 = -inv(A11)*A12*inv(A22)
 */
#define qmod(a,b) ((a)-(__mul24((b),(a)/(b))))
__global__ void
triple_dgemm_update_16_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
//	const int page = (blockIdx.y)%(npages);
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * (blockDim.x*blockDim.y);
	const int iby = bIdy * 16;
	const int id = inx + iny*blockDim.x;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part one---------------------------//
	{
		// A12*inv(A22) -> A12
		// A=A12, B=inv(A22), C=A12(d_dinvA)
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
				   (qmod(page,PagesPerNB))*(blk*2)*NB +
				   (qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + blk*lda + page*blk*2;  
		B = d_dinvA + blk*NB + blk;
		C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+4][iny]    = B[4+0*ldb];
			bs[inx+4][iny+4]  = B[4+4*ldb];
			bs[inx+4][iny+8]  = B[4+8*ldb];
			bs[inx+4][iny+12] = B[4+12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			bs[inx+12][iny]    = B[12+0*ldb];
			bs[inx+12][iny+4]  = B[12+4*ldb];
			bs[inx+12][iny+8]  = B[12+8*ldb];
			bs[inx+12][iny+12] = B[12+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
			
	__syncthreads();
	
	//--------------------------part two---------------------------//
	{
		// -inv(A11)*A12 -> A12
		// A=inv(A11), B=A12, C=A12
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		A = d_dinvA;
		B = C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+4][iny]    = B[4+0*ldb];
			bs[inx+4][iny+4]  = B[4+4*ldb];
			bs[inx+4][iny+8]  = B[4+8*ldb];
			bs[inx+4][iny+12] = B[4+12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			bs[inx+12][iny]    = B[12+0*ldb];
			bs[inx+12][iny+4]  = B[12+4*ldb];
			bs[inx+12][iny+8]  = B[12+8*ldb];
			bs[inx+12][iny+12] = B[12+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
}

/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
#define qmod(a,b) ((a)-(__mul24((b),(a)/(b))))
__global__ void
triple_dgemm_update_16_part1_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
//	const int page = (blockIdx.y)%(npages);
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * (blockDim.x*blockDim.y);
	const int iby = bIdy * 16;
	const int id = inx + iny*blockDim.x;
	__shared__ double bs[16][17];

	//--------------------------part one---------------------------//
	{
		// A21*inv(A11) -> A21
		// A=A21, B=inv(A11), C=A21
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		int PagesPerNB = NB/(blk*2);

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + page*blk*2 + blk;  
		B = d_dinvA; 
		C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+4][iny]    = B[4+0*ldb];
			bs[inx+4][iny+4]  = B[4+4*ldb];
			bs[inx+4][iny+8]  = B[4+8*ldb];
			bs[inx+4][iny+12] = B[4+12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			bs[inx+12][iny]    = B[12+0*ldb];
			bs[inx+12][iny+4]  = B[12+4*ldb];
			bs[inx+12][iny+8]  = B[12+8*ldb];
			bs[inx+12][iny+12] = B[12+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
			
	__syncthreads();
	
}

/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
#define qmod(a,b) ((a)-(__mul24((b),(a)/(b))))
__global__ void
triple_dgemm_update_16_part2_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * (blockDim.x*blockDim.y);
	const int iby = bIdy * 16;
	const int id = inx + iny*blockDim.x;
	__shared__ double bs[16][17];
	
	//--------------------------part two---------------------------//
	{
		// -inv(A22)*A21 -> A21
		// A=inv(A22), B=A21, C=A21
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		int PagesPerNB = NB/(blk*2);
		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = d_dinvA + blk*NB + blk;
		B = C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+4][iny]    = B[4+0*ldb];
			bs[inx+4][iny+4]  = B[4+4*ldb];
			bs[inx+4][iny+8]  = B[4+8*ldb];
			bs[inx+4][iny+12] = B[4+12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			bs[inx+12][iny]    = B[12+0*ldb];
			bs[inx+12][iny+4]  = B[12+4*ldb];
			bs[inx+12][iny+8]  = B[12+8*ldb];
			bs[inx+12][iny+12] = B[12+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
	__syncthreads();
}

/*
 * B21 = -inv(A11)*A12*inv(A22)
 */
__global__ void
triple_dgemm_update_32_part1_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * (blockDim.x*blockDim.y);
	const int iby = bIdy * 16;
	const int id = inx + iny*blockDim.x;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part one---------------------------//
	{
		// A12*inv(A22) -> A21
		// A=A12, B=inv(A22), C=A12(d_dinvA)
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
				   (qmod(page,PagesPerNB))*(blk*2)*NB +
				   (qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + blk*lda + page*blk*2;  
		B = d_dinvA + blk*NB + blk;
		C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
			
	__syncthreads();
}
/*
 * B21 = -inv(A11)*A12*inv(A22)
 */
__global__ void
triple_dgemm_update_32_part2_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * (blockDim.x*blockDim.y);
	const int iby = bIdy * 16;
	const int id = inx + iny*blockDim.x;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	
	//--------------------------part two---------------------------//
	{
		// -inv(A11)*A12 -> A12
		// A=inv(A11), B=A12, C=A12
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
				   (qmod(page,PagesPerNB))*(blk*2)*NB +
				   (qmod(page,PagesPerNB))*(blk*2);

		A = d_dinvA;
		B = C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
}
/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
__global__ void
triple_dgemm_update_32_part1_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * (blockDim.x*blockDim.y);
	const int iby = bIdy * 16;
	const int id = inx + iny*blockDim.x;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part one---------------------------//
	{
		// A21*inv(A11) -> A21
		// A=A21, B=inv(A11), C=A21
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + page*blk*2 + blk;  
		B = d_dinvA;
		C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
			
	__syncthreads();
}
/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
__global__ void
triple_dgemm_update_32_part2_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * (blockDim.x*blockDim.y);
	const int iby = bIdy * 16;
	const int id = inx + iny*blockDim.x;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part two---------------------------//
	{
		// -inv(A22)*A21 -> A21
		// A=inv(A22), B=A21, C=A21
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = d_dinvA + blk*NB + blk;
		B = C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			bs[inx+8][iny]    = B[8+0*ldb];
			bs[inx+8][iny+4]  = B[8+4*ldb];
			bs[inx+8][iny+8]  = B[8+8*ldb];
			bs[inx+8][iny+12] = B[8+12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
}

/*
 * B21 = -inv(A11)*A12*inv(A22)
 */
__global__ void
triple_dgemm_update_64_part1_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part one---------------------------//
	{
		// A12*inv(A22) -> A12(d_dinvA)
		// A=A12, B=inv(A22), C=A12
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
				   (qmod(page,PagesPerNB))*(blk*2)*NB +
				   (qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + blk*lda + page*blk*2;  
		B = d_dinvA + blk*NB + blk;
		C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
			
}
/*
 * B21 = -inv(A11)*A12*inv(A22)
 */
__global__ void
triple_dgemm_update_64_part2_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
			
	//--------------------------part two---------------------------//
	{
		// -inv(A11)*A12 -> A12
		// A=inv(A11), B=A12, C=A12
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
				   (qmod(page,PagesPerNB))*(blk*2)*NB +
				   (qmod(page,PagesPerNB))*(blk*2);

		A = d_dinvA;
		B = C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
}
/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
__global__ void
triple_dgemm_update_64_part1_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part one---------------------------//
	{
		// A21*inv(A11) -> A21
		// A=A21, B=inv(A11), C=A21
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + page*blk*2 + blk;  
		B = d_dinvA;
		C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
}
/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
__global__ void
triple_dgemm_update_64_part2_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
			
	//--------------------------part two---------------------------//
	{
		// -inv(A22)*A21 -> A21
		// A=inv(A22), B=A21, C=A21
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = d_dinvA + blk*NB + blk;
		B = C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
}

/*
 * B21 = -inv(A11)*A12*inv(A22)
 */
__global__ void
triple_dgemm_update_above64_part1_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part one---------------------------//
	{
		// A12*inv(A22) -> A12(d_dinvA)
		// A=A12, B=inv(A22), C=A12
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
				   (qmod(page,PagesPerNB))*(blk*2)*NB +
				   (qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + blk*lda + page*blk*2;  
		B = d_dinvA + blk*NB + blk;
		C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
			
}
/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
__global__ void
triple_dgemm_update_above64_part1_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	//--------------------------part one---------------------------//
	{
		// A21*inv(A11) -> A21
		// A=A21, B=inv(A11), C=A21
		double *A, *B, *C;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = Ain + page*lda*blk*2 + page*blk*2 + blk;  
		B = d_dinvA;
		C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = c[i];
	}
			
}

/*
 * B21 = -inv(A11)*A12*inv(A22)
 */
__global__ void
triple_dgemm_update_above64_part2_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	
	//--------------------------part two---------------------------//
	{
		// -inv(A11)*A12 -> A12
		// A=inv(A11), B=A12, C=A12
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
				   (qmod(page,PagesPerNB))*(blk*2)*NB +
				   (qmod(page,PagesPerNB))*(blk*2);

		A = d_dinvA;
		B = d_dinvA + blk*NB;
		C = d_dinvA + blk;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
}
/*
 * part 3, copy data into position 
 */
__global__ void
triple_dgemm_update_above64_part3_R (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;

	int PagesPerNB = NB/(blk*2);
	
	//--------------------------part two---------------------------//
	{
		// -inv(A11)*A12 -> A12
		// A=inv(A11), B=A12, C=A12
		double *C_temp, *C_real;
		int ldc = NB;

		C_temp = d_dinvA + NB*NB*(page/PagesPerNB) + 
					  (qmod(page,PagesPerNB))*(blk*2)*NB +
					  (qmod(page,PagesPerNB))*(blk*2) +
					  blk;
		C_real = d_dinvA + NB*NB*(page/PagesPerNB) + 
					  (qmod(page,PagesPerNB))*(blk*2)*NB +
					  blk*NB +
					  (qmod(page,PagesPerNB))*(blk*2);

		C_temp += ibx + id  + __mul24( iby, ldc );
		C_real += ibx + id  + __mul24( iby, ldc );


		for( int i = 0; i < 16; i++, C_temp+=ldc, C_real+=ldc )
		{
			C_real[0] = C_temp[0];
			C_temp[0] = 0;
		}

	}
}

/*
 * part 3: copy data back to position 
 */
__global__ void
triple_dgemm_update_above64_part3_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;

	int PagesPerNB = NB/(blk*2);
	
	//--------------------------part three---------------------------//
	{
		// -inv(A22)*A21 -> A21
		// A=inv(A22), B=A21, C=A21
		double *C_temp, *C_real;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		C_real = d_dinvA + blk;
		
		C_temp = d_dinvA + blk*NB;

		C_temp += ibx + id  + __mul24( iby, ldc );
		C_real += ibx + id  + __mul24( iby, ldc );

		for( int i = 0; i < 16; i++, C_real+=ldc, C_temp+=ldc )
		{
			C_real[0] = C_temp[0];
			C_temp[0] = 0;
		}
	}
	__syncthreads();
}
/*
 * B21 = -inv(A22)*A21*inv(A11)
 */
__global__ void
triple_dgemm_update_above64_part2_L (double * Ain, double *d_dinvA, int blk, int lda, int npages)
{
	const int bIdy = blockIdx.y/npages;
	const int page = qmod(blockIdx.y, npages);
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x*64;
	const int iby = bIdy*16;
	const int id = inx + iny*16;
	__shared__ double bs[16][17];

	int PagesPerNB = NB/(blk*2);
	
	//--------------------------part two---------------------------//
	{
		// -inv(A22)*A21 -> A21
		// A=inv(A22), B=A21, C=A21
		double *A, *B, *C;
		int lda = NB;
		int ldb = NB;
		int ldc = NB;

		d_dinvA += NB*NB*(page/PagesPerNB) + 
					(qmod(page,PagesPerNB))*(blk*2)*NB +
					(qmod(page,PagesPerNB))*(blk*2);

		A = d_dinvA + blk*NB + blk;
		B = d_dinvA + blk;
		
		C = d_dinvA + blk*NB;

		A += ibx + id;
		B += inx + __mul24( iby + iny, ldb );
		C += ibx + id  + __mul24( iby, ldc );

		const double *Blast = B + blk;

		double c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		do
		{
			double a[4] = { A[0*lda], A[1*lda], A[2*lda], A[3*lda] };

			bs[inx][iny]    = B[0*ldb];
			bs[inx][iny+4]  = B[4*ldb];
			bs[inx][iny+8]  = B[8*ldb];
			bs[inx][iny+12] = B[12*ldb];
			__syncthreads();

			A += 4*lda;
			ssaxpy( a[0], &bs[0][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[1][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[2][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[3][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[4][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[5][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[6][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[7][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[8][0], c );		a[0] = A[0*lda];
			ssaxpy( a[1], &bs[9][0], c );		a[1] = A[1*lda];
			ssaxpy( a[2], &bs[10][0], c );		a[2] = A[2*lda];
			ssaxpy( a[3], &bs[11][0], c );		a[3] = A[3*lda];

			A += 4*lda;
			ssaxpy( a[0], &bs[12][0], c );
			ssaxpy( a[1], &bs[13][0], c );
			ssaxpy( a[2], &bs[14][0], c );
			ssaxpy( a[3], &bs[15][0], c );

			B += 16;
			__syncthreads();
		} while( B < Blast );

		for( int i = 0; i < 16; i++, C += ldc )
			C[0] = (-1)*c[i];
	}
}

__global__ void
b_copy_kernel (int M, int N, double *b, int ldb, double *d_x, int ldx)
{
	int by = blockIdx.y;

	int gx = blockIdx.x*blockDim.x+threadIdx.x;

	if (gx < M)
		b[by*ldb+gx] = d_x[by*ldx+gx];
}


extern "C"
void diag_dtrtri (int M, char uplo, char diag, double *A, double *d_dinvA, int lda)
{
	int nblocks = M/BLOCK_SIZE+(M%BLOCK_SIZE!=0);

	if (uplo == 'l' || uplo == 'L')
	{
		// solve the diagonal blocks
		diag_dtrtri_kernel_lower<<<nblocks, BLOCK_SIZE>>>(diag, A, d_dinvA, lda);

		// update the inverse up to the size of BLOCK_SIZE
		for (int i=BLOCK_SIZE; i<NB; i*=2)
		{
			int npages = M/(i*2)+(M%(i*2)!=0);
			dim3 dimBlock((i<=32)?(i/4):16, 4);
			dim3 dimGrid(i/(dimBlock.x*dimBlock.y), npages*(i/16));	//emulated 3D grid, see 3d_grid.txt 
			
			switch (i)
			{
				case 16:
					triple_dgemm_update_16_part1_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_16_part2_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
				case 32:
					triple_dgemm_update_32_part1_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_32_part2_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
				case 64:
					triple_dgemm_update_64_part1_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_64_part2_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
				default:
					triple_dgemm_update_above64_part1_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_above64_part2_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_above64_part3_L<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
			}
			if (i*2>=M) break;
		}
	}
	else
	{
		diag_dtrtri_kernel_upper<<<nblocks, BLOCK_SIZE>>>(diag, A, d_dinvA, lda);

		// update the inverse up to the size of BLOCK_SIZE
		for (int i=BLOCK_SIZE; i<NB; i*=2)
		{
			int npages = M/(i*2)+(M%(i*2)!=0);
			dim3 dimBlock((i<=32)?(i/4):16, 4);
			dim3 dimGrid(i/(dimBlock.x*dimBlock.y), npages*(i/16));	//emulated 3D grid, see 3d_grid.txt 
			
			switch (i)
			{
				case 16:
					triple_dgemm_update_16_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
				case 32:
					triple_dgemm_update_32_part1_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_32_part2_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
				case 64:
					triple_dgemm_update_64_part1_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_64_part2_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
				default:
					triple_dgemm_update_above64_part1_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_above64_part2_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					triple_dgemm_update_above64_part3_R<<<dimGrid, dimBlock>>>(A, d_dinvA, i, lda, npages);
					break;
			}
			if (i*2>=M) break;
		}

	}
}

/*
 * magmablas_dtrsm
 */
extern "C"
void magmablas_dtrsm_tesla( char side, char uplo, char tran, char diag, int M, int N, 
                            double alpha, /*const*/ double* A, int lda, double* b, int ldb)
{
	/*  -- magma (version 0.1) --
		univ. of tennessee, knoxville
		univ. of california, berkeley
		univ. of colorado, denver
		october 2009

		purpose
		=======

		dtrsm  solves one of the matrix equations on gpu

		op( a )*x = alpha*b,   or   x*op( a ) = alpha*b,

		where alpha is a scalar, x and b are m by n matrices, a is a unit, or
		non-unit,  upper or lower triangular matrix  and  op( a )  is one  of

		op( A ) = A   or   op( A ) = A'.

		The matrix X is overwritten on B.

		When M or N is not a multiple of blocking size, which is 32 for now, cublasDtrsm will
		be called instead. There soon will not be this limitation both for arbitrary problem 
		size and blocking size.
	   
	   
	   Arguments
	   ==========
	   
	   side   - CHARACTER*1.
	            On entry, side specifies whether op( A ) appears on the left
	            or right of X as follows:
	   
	               side = 'L' or 'l'   op( A )*X = alpha*B.
	   
	               side = 'R' or 'r'   X*op( A ) = alpha*B.
	   
	            Unchanged on exit.
	   
	   uplo   - CHARACTER*1.
	            On entry, uplo specifies whether the matrix A is an upper or
	            lower triangular matrix as follows:
	   
	               uplo = 'U' or 'u'   A is an upper triangular matrix.
	   
	               uplo = 'L' or 'l'   A is a lower triangular matrix.
	   
	            Unchanged on exit.
	   
	   tran - CHARACTER*1.
	            On entry, tran specifies the form of op( A ) to be used in
	            the matrix multiplication as follows:
	   
	               tran = 'N' or 'n'   op( A ) = A.
	   
	               tran = 'T' or 't'   op( A ) = A'.
	   
	               tran = 'C' or 'c'   op( A ) = A'.
	   
	            Unchanged on exit.
	   
	   diag   - CHARACTER*1.
	            On entry, diag specifies whether or not A is unit triangular
	            as follows:
	   
	               diag = 'U' or 'u'   A is assumed to be unit triangular.
	   
	               diag = 'N' or 'n'   A is not assumed to be unit
	                                   triangular.
	   
	            Unchanged on exit.
	   
	   m      - INTEGER.
	            On entry, m specifies the number of rows of B. m must be at
	            least zero.
	            Unchanged on exit.
	   
				n      - INTEGER.
	             On entry, n specifies the number of columns of B.  n must be
	             at least zero.
	             Unchanged on exit.
	   
	    alpha  - REAL.
	             On entry,  alpha specifies the scalar  alpha. When  alpha is
	             zero then  A is not referenced and  B need not be set before
	             entry.
	             Unchanged on exit.
	   
	    A      - REAL             array of DIMENSION ( lda, k ), where k is m
	             when  side = 'L' or 'l'  and is  n  when  side = 'R' or 'r'.
	             Before entry  with  uplo = 'U' or 'u',  the  leading  k by k
	             upper triangular part of the array  A must contain the upper
	             triangular matrix  and the strictly lower triangular part of
	             A is not referenced.
	             Before entry  with  uplo = 'L' or 'l',  the  leading  k by k
	             lower triangular part of the array  A must contain the lower
	             triangular matrix  and the strictly upper triangular part of
	             A is not referenced.
	             Note that when  diag = 'U' or 'u',  the diagonal elements of
	             A  are not referenced either,  but are assumed to be  unity.
	             Unchanged on exit.
	   
	    lda    - INTEGER.
	             On entry, lda specifies the first dimension of A as declared
	             in the calling (sub) program.  When  side = 'L' or 'l'  then
	             lda  must be at least  max( 1, m ),  when  side = 'R' or 'r'
	             then lda must be at least max( 1, n ).
	             Unchanged on exit.
	   
	    b      - REAL             array of DIMENSION ( ldb, n ).
	             Before entry,  the leading  m by n part of the array  B must
	             contain  the  right-hand  side  matrix  B,  and  on exit  is
	             overwritten by the solution matrix  X.
	   
	    ldb    - INTEGER.
	             On entry, ldb specifies the first dimension of B as declared
	             in  the  calling  (sub)  program.   ldb  must  be  at  least
	             max( 1, m ).
	             Unchanged on exit.
	   
	   
	    Level 3 Blas routine.
		*
    ===================================================================== */

	int i;
	double *d_dinvA, *d_x;

	/* quick return on wrong size */
	if (M<=0 || N<=0)
		return;

	if (side == 'l' || side == 'L')
	{
		/* inverse the diagonals
		 * Allocate device memory for the inversed diagonal blocks, size=m*NB
		 */
		cudaMalloc((void**)&d_dinvA, NB*((M/NB)+(M%NB!=0))*NB*sizeof(double));
		cudaMalloc((void**)&d_x, N*M*sizeof(double));
		cudaMemset(d_x, 0, N*M*sizeof(double));
		cudaMemset (d_dinvA, 0, NB*((M/NB)+(M%NB!=0))*NB*sizeof(double));
		diag_dtrtri (M, uplo, diag, A, d_dinvA, lda);

		if (tran == 'N' || tran == 'n')
		/* the non-transpose case */
		{
			if (uplo == 'L' || uplo == 'l')
			{
			/* the lower case */
				
				/* handle the first block seperately with alpha */

				int MM = min (NB, M); 
				cublasDgemm ('N', 'N', MM, N, MM, alpha, d_dinvA, NB, b, ldb, 0, d_x, M);  

				if (NB>=M)
				{
					b_copy();
					cudaFree(d_dinvA);
					cudaFree(d_x);
					return;
				}

				cublasDgemm ('N', 'N', M-NB, N, NB, -1.0, A+NB, lda, d_x, M, alpha, b+NB, ldb);

				/* the rest blocks */
				for (i=NB; i<M; i+=NB)
				{
					MM = min (M-i, NB);
					cublasDgemm ('N', 'N', MM, N, MM, 1.0, d_dinvA+i*NB, NB, b+i, ldb, 0, d_x+i, M);  
					
					if (i+NB>=M)
						break;

					cublasDgemm ('N', 'N', M-i-NB, N, NB, -1.0, A+i*lda+i+NB, lda, d_x+i, M, 1.0, b+i+NB, ldb);
				}
			}
			else
			{
			/* the upper case */

				/* handle the first block seperately with alpha */
				int MM = (M%NB==0)?NB:(M%NB); 
				i = M-MM;
				cublasDgemm ('N', 'N', MM, N, MM, alpha, d_dinvA+i*NB, NB, b+i, ldb, 0.0, d_x+i, M); 
					
				if (i-NB<0)
				{
					b_copy();
					cudaFree(d_dinvA);
					cudaFree(d_x);
					return;
				}

				cublasDgemm ('N', 'N', i, N, MM, -1.0, A+i*lda, lda, d_x+i, M, alpha, b, ldb);

				/* the rest blocks */
				for (i=M-MM-NB; i>=0; i-=NB)
				{
					cublasDgemm ('N', 'N', NB, N, NB, 1.0, d_dinvA+i*NB, NB, b+i, ldb, 0.0, d_x+i, M);

					if (i-NB<0)
						break;

					cublasDgemm ('N', 'N', i, N, NB, -1.0, A+i*lda, lda, d_x+i, M, 1.0, b, ldb);
				}

			}
		}
		else
		/* the transpose case */
		{
			if (uplo == 'L' || uplo == 'l')
			{
			/* the lower case */
				
				/* handle the first block seperately with alpha */
				int MM = (M%NB==0)?NB:(M%NB); 
				i=M-MM; 
				cublasDgemm ('T', 'N', MM, N, MM, alpha, d_dinvA+i*NB, NB, b+i, ldb, 0, d_x+i, M);  

				if (i-NB<0)
				{
					b_copy();
					cudaFree(d_dinvA);
					cudaFree(d_x);
					return;
				}

				cublasDgemm ('T', 'N', i, N, MM, -1.0, A+i, lda, d_x+i, M, alpha, b, ldb);

				/* the rest blocks */
				for (i=M-MM-NB; i>=0; i-=NB)
				{
					cublasDgemm ('T', 'N', NB, N, NB, 1.0, d_dinvA+i*NB, NB, b+i, ldb, 0, d_x+i, M);  

					if (i-NB<0)
						break;

					cublasDgemm ('T', 'N', i, N, NB, -1.0, A+i, lda, d_x+i, M, 1.0, b, ldb);
				}
			}
			else
			{
			/* the upper case */
					
				/* handle the first block seperately with alpha */
				int MM = min (NB, M); 
				cublasDgemm ('T', 'N', MM, N, MM, alpha, d_dinvA, NB, b, ldb, 0, d_x, M);  

				if (NB>=M)
				{
					b_copy();
					cudaFree(d_dinvA);
					cudaFree(d_x);
					return;
				}

				cublasDgemm ('T', 'N', M-NB, N, NB, -1.0, A+(NB)*lda, lda, d_x, M, alpha, b+NB, ldb);

				/* the rest blocks */
				for (i=NB; i<M; i+=NB)
				{
					MM = min (M-i, NB);
					cublasDgemm ('T', 'N', MM, N, MM, 1.0, d_dinvA+i*NB, NB, b+i, ldb, 0, d_x+i, M);  
					
					if (i+NB>=M)
						break;

					cublasDgemm ('T', 'N', M-i-NB, N, NB, -1.0, A+(i+NB)*lda+i, lda, d_x+i, M, 1.0, b+i+NB, ldb);
				}
			}
		}
	}
	else
	{	// side=R

		/* inverse the diagonals
		 * Allocate device memory for the inversed diagonal blocks, size=N*BLOCK_SIZE 
		 */
		cudaMalloc((void**)&d_dinvA, NB*((N/NB)+(N%NB!=0))*NB*sizeof(double));
		cudaMalloc((void**)&d_x, N*M*sizeof(double));
		cudaMemset(d_x, 0, N*M*sizeof(double));
		cudaMemset (d_dinvA, 0, NB*((N/NB)+(N%NB!=0))*NB*sizeof(double));
		diag_dtrtri (N, uplo, diag, A, d_dinvA, lda);

		if (tran == 'N' || tran == 'n')
		/* the non-transpose case */
		{
			if (uplo == 'L' || uplo == 'l')
			{
			/* the lower case */
				
				/* handle the first block seperately with alpha */
				int NN = (N%NB==0)?NB:(N%NB);
				i=N-NN;
				cublasDgemm ('N', 'N', M, NN, NN, alpha, b+ldb*i, ldb, d_dinvA+i*NB, NB, 0.0, d_x+i*M, M); 

				if (i-NB<0)
				{
					b_copy();
					cudaFree(d_x);
					cudaFree(d_dinvA);
					return;
				}

				cublasDgemm ('N', 'N', M, i, NN, -1.0, d_x+i*M, M, A+i, lda, alpha, b, ldb);

				/* the rest blocks */
				for (i=N-NN-NB; i>=0; i-=NB)
				{
					cublasDgemm ('N', 'N', M, NB, NB, 1.0, b+ldb*i, ldb, d_dinvA+i*NB, NB, 0.0, d_x+i*M, M); 
					
					if (i-NB<0)
						break;

					cublasDgemm ('N', 'N', M, i, NB, -1.0, d_x+i*M, M, A+i, lda, 1.0, b, ldb);
				}
			}
			else
			{
			/* the upper case */
				
				/* handle the first block seperately with alpha */
				int NN = min(NB, N); 
				cublasDgemm ('N', 'N', M, NN, NN, alpha, b, ldb, d_dinvA, NB, 0, d_x, M);  

				if (NB>=N)
				{
					b_copy();
					cudaFree(d_x);
					cudaFree(d_dinvA);
					return;
				}

				cublasDgemm ('N', 'N', M, N-NB, NB, -1.0, d_x, M, A+NB*lda, lda, alpha, b+NB*ldb, ldb);
				
				/* the rest blocks */
				for (i=NB; i<N; i+=NB)
				{
					NN = min(NB, N-i); 
					cublasDgemm ('N', 'N', M, NN, NN, 1.0, b+ldb*i, ldb, d_dinvA+i*NB, NB, 0, d_x+i*M, M);  

					if (i+NB>=N)
						break;

					cublasDgemm ('N', 'N', M, N-i-NB, NB, -1.0, d_x+i*M, M,   A+(i+NB)*lda+i, lda, 1.0, b+(i+NB)*ldb, ldb);
				}
			}
		}
		else
		/* the transpose case */
		{
			if (uplo == 'L' || uplo == 'l')
			{
			/* the lower case */
				
				/* handle the first block seperately with alpha */
				int NN = min(NB, N); 
				cublasDgemm ('N', 'T', M, NN, NN, alpha, b, ldb, d_dinvA, NB, 0, d_x, M);  

				if (NB>=N)
				{
					b_copy();
					cudaFree(d_x);
					cudaFree(d_dinvA);
					return;
				}

				cublasDgemm ('N', 'T', M, N-NB, NB, -1.0, d_x, M, A+NB, lda, alpha, b+NB*ldb, ldb);

				/* the rest blocks */
				for (i=NB; i<N; i+=NB)
				{
					NN = min(NB, N-i); 
					cublasDgemm ('N', 'T', M, NN, NN, 1.0, b+ldb*i, ldb, d_dinvA+i*NB, NB, 0, d_x+i*M, M);  

					if (i+NB>=N)
						break;

					cublasDgemm ('N', 'T', M, N-i-NB, NB, -1.0, d_x+i*M, M,   A+i*lda+NB+i, lda, 1.0, b+(i+NB)*ldb, ldb);
				}
			}
			else
			{
			/* the upper case */
				
				/* handle the first block seperately with alpha */
				int NN = (N%NB==0)?NB:(N%NB);
				i=N-NN;
				cublasDgemm ('N', 'T', M, NN, NN, alpha, b+ldb*i, ldb, d_dinvA+i*NB, NB, 0.0, d_x+i*M, M); 

				if (i-NB<0)
				{
					b_copy();
					cudaFree(d_x);
					cudaFree(d_dinvA);
					return;
				}

				cublasDgemm ('N', 'T', M, i, NN, -1.0, d_x+i*M, M, A+i*lda, lda, alpha, b, ldb);
				
				/* the rest blocks */
				for (i=N-NN-NB; i>=0; i-=NB)
				{
					cublasDgemm ('N', 'T', M, NB, NB, 1.0, b+ldb*i, ldb, d_dinvA+i*NB, NB, 0.0, d_x+i*M, M); 

					if (i-NB<0)
						break;

					cublasDgemm ('N', 'T', M, i, NB, -1.0, d_x+i*M, M, A+i*lda, lda, 1.0, b, ldb);
				}
			}
		}
	}
		
	b_copy();
	cudaFree(d_dinvA);
	cudaFree(d_x);
}

