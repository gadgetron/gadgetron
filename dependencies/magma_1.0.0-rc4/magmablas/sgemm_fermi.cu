/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/


/*
    blk_M=96 blk_N=96 blk_K=16 nthd_x=64 nthd_y=4
*/
#include "common_magma.h"
#include "commonblas_s.h"

#define magmablas_sgemm_fermi magmablas_sgemm

texture<float,1>  tex_x_float_A;
texture<float,1>  tex_x_float_B;

static __inline__ __device__ float fetch_x_A(const int& i, const float * x)
{
	return tex1Dfetch(tex_x_float_A, i);
}

static __inline__ __device__ float fetch_x_B(const int& i, const float * x)
{
	return tex1Dfetch(tex_x_float_B, i);
}

extern "C" __global__ void 
fermiSgemm_v2_kernel_NN(float *C, const float *A, const float *B,  
                        int m, int n,  int k,  int lda,  int ldb,  int ldc, 
                        float alpha, float beta,
                        int offsetA, int offsetB) 
{
	const  int tx = threadIdx.x;
	const  int ty = threadIdx.y;

	const int iby = blockIdx.y * 96;
	const int ibx = blockIdx.x * 96;
	const int idt = ty * 64 + tx;

	const int tx2= idt%16;	// idx2
	const int ty2= idt/16;	// idy2

	__shared__ float Bb[16][97];
	__shared__ float Abs[96][17];

	float xxA[6];
	float xxB[6];
	
	int trackA = offsetA + ibx + tx2 + __mul24(ty2, lda);
	int trackB = offsetB + tx2 + __mul24(iby+ty2*6, ldb);

	int tll = ty2;
	A += trackA; 
	B += trackB; 

	// read a block of 96x16 to A and 16x96 to B
	// each thread reads 6 data point, 1 point in each 16x16 subblock
	#pragma unroll
	for(int y=0; y<6; y++)
		Abs[tx2+y*16][ty2] = /* (tll<k)* */ fetch_x_A(trackA + y*16, A);

	#pragma unroll
	for(int y=0; y<6; y++)
		Bb[tx2][ty2*6+y] = fetch_x_B( trackB + y*ldb, B) ;

	__syncthreads();

	const float *Bend = B + k-16;

	float Axs[6];
	float Bxp[6];

	float Cb[36] = {0,0,0,0,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, 0,0,0,0,0,0,
		0,0,0,0,0,0, 0,0,0,0,0,0};

	do 
	{
		tll += 16;
		A += lda*16;
		B += 16;
		trackA+=16*lda ; 
		trackB+=16;

		// calculate part of C using the first 96x8 of A and 8x96 of B
		#pragma unroll 
		for( int j1=0;j1<8;j1++)
		{
			#pragma unroll
			for( int y=0; y<6; y++)
				Axs[y] = Abs[tx2+y*16][j1];

			#pragma unroll
			for( int y=0; y<6; y++)
				Bxp[y]= Bb[j1][ty2+y*16];

			#pragma unroll 
			for( int x=0; x<6; x++)
			{
				#pragma unroll 
				for( int y=0; y<6; y++)
				{
					Cb[x*6 + y] += Axs[x]*Bxp[y];
				}
			}
		}

		// pre-read the next A and B
		#pragma unroll
		for( int y=0; y<6; y++)
			xxA[y] = /* (tll<k)* */ fetch_x_A(trackA + y*16, A);	
			// without going through texture and tll protection, 
			// nonzeros are fetched
			// texture boundary control seems to be providing
			// safety here but officially tex1Dfetch is not suppoted
			// by clamping or filtering (programming guide B.8.1)

		#pragma unroll
		for( int y=0; y<6; y++)
			xxB[y] = fetch_x_B( trackB + y*ldb, B);

		// calculate another part of C using the 2nd 96x8 of A and 8x96 of B
		#pragma unroll 
		for( int j1=8;j1<16;j1++)
		{
			#pragma unroll
			for( int y=0;y<6;y++)
				Axs[y] = Abs[tx2+y*16][j1] ;

			#pragma unroll
			for( int y=0;y<6;y++)
				Bxp[y]= Bb[j1][ty2+y*16];

			#pragma unroll 
			for( int x=0;x<6;x++)
			{
				#pragma unroll 
				for( int y=0; y<6; y++)
				{
					Cb[x*6+y] += Axs[x]*Bxp[y];
				}
			}
		}

		__syncthreads();

		// put the next A and B into position
		#pragma unroll
		for(int y=0;y<6;y++)
			Abs[tx2+y*16][ty2] =xxA[y];

		#pragma unroll
		for(int y=0; y<6; y++)
			Bb[tx2][ty2*6 + y] =xxB[y];

		__syncthreads();
	} 
	while (B < Bend);

	// C += ty2 + ibx  + __mul24 (tx2 +  iby ,ldc);
	C += tx2 + ibx  + __mul24 (ty2 +  iby ,ldc);

	// tail case
	#pragma unroll 
	for( int j1=0;j1<16;j1++)
	{

		#pragma unroll
		for( int y=0;y<6;y++)
			Bxp[y]= Bb[j1][ty2+y*16];

		#pragma unroll
		for( int y=0;y<6;y++)
			Axs[y] = Abs[tx2+y*16][j1] ;

		#pragma unroll 
		for( int x=0;x<6;x++)
			#pragma unroll 
			for( int y=0; y<6; y++)
				Cb[x*6+y]  += Axs[x]*Bxp[y];	
	}

	// __syncthreads();

	// C += ty2 + ibx  + __mul24 (tx2 +  iby ,ldc);
	int gy = iby + ty2;
	// writing C
	#pragma unroll
	for( int y=0; y<6; y++, gy+=16)
	{
		int gx = ibx + tx2; 
		#pragma unroll
		for(int x=0; x<6; x++, gx+=16)
		{
			if (gx < m && gy < n)
				C[x*16] = alpha*Cb[y+x*6] + beta * C[x*16];
		}

		C+=ldc*16;
	}
}

//========================================================================

extern "C" __global__ void 
fermiSgemm_v2_kernel_TN(float *C, const float *A, const float *B,  
                        int m, int n,  int k,  int lda,  int ldb,  
                        int ldc, float alpha, float beta,
                        int offsetA, int offsetB) 
{
	const  int tx = threadIdx.x;
	const  int ty = threadIdx.y;

	const int iby = blockIdx.y * 96;
	const int ibx = blockIdx.x * 96;
	const int idt = ty * 64 + tx;

	const int tx2 = idt%16;
	const int ty2 = idt/16;

	__shared__ float Bb[16][97];
	__shared__ float Abs[96][17];

	float xxA[6];
	float xxB[6];
	
	int trackA = offsetA + tx2+ __mul24(ibx + ty2*6, lda );
	int trackB = offsetB + tx2+ __mul24(iby + ty2*6, ldb );
	
	A+= trackA; 
	B+= trackB; 

	int tll = tx2;

	#pragma unroll
	for(int y=0; y<6; y++)
		Abs[ty2*6+y][tx2] = (tll<k)*fetch_x_A(trackA + y*lda, A);

	#pragma unroll
	for(int y=0; y<6; y++)
		Bb[tx2][ty2*6+y] = /* (tll<k)* */ fetch_x_B(trackB + y*ldb, B ) ;

	__syncthreads();

	const float *Bend = B + k-16;

	float Axs[6];
	float Bxp[6];

	float Cb[36] = {0,0,0,0,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, 0,0,0,0,0,0,
		0,0,0,0,0,0, 0,0,0,0,0,0};
	do 
	{
		tll += 16;
		A += 16;
		B += 16;
		trackA += 16; 
		trackB += 16;

		// pre-read the next strip of A and B
		#pragma unroll
		for( int y=0; y<6; y++)
			xxA[y] = (tll<k)*fetch_x_A(trackA + y*lda, A);

		#pragma unroll
		for( int y=0; y<6; y++)
			xxB[y] = /* (tll<k)* */ fetch_x_B(trackB + y*ldb, B);

		// computing
		#pragma unroll 
		for( int j1=0; j1<16; j1++)
		{
			#pragma unroll
			for( int y=0; y<6; y++)
				Axs[y] = Abs[tx2 + y*16][j1];

			#pragma unroll
			for( int y=0; y<6; y++)
				Bxp[y]= Bb[j1][ty2 + y*16];

			#pragma unroll 
			for( int x=0;x<6; x++)
				#pragma unroll 
				for(int y=0; y<6; y++)
					Cb[x*6 + y] += Axs[x]*Bxp[y];
		}
		__syncthreads();
		
		#pragma unroll
		for(int y=0; y<6; y++)
			Abs[ty2*6 + y][tx2] = xxA[y]; 

		#pragma unroll
		for( int y=0; y<6; y++)
			Bb[tx2][ty2*6 + y] = xxB[y];
		__syncthreads();
	} 
	while (B < Bend);

	C += tx2 + ibx  + __mul24 (ty2 + iby, ldc);

	#pragma unroll 
	for( int j1=0; j1<16; j1++)
	{

		#pragma unroll
		for( int y=0; y<6; y++)
			Axs[y] =  Abs[tx2 + y*16][j1] ;

		#pragma unroll
		for( int y=0; y<6; y++)
			Bxp[y] = Bb[j1][ty2 + y*16];

		#pragma unroll 
		for( int x=0; x<6; x++)
			#pragma unroll 
			for( int y=0; y<6; y++)
				Cb[x*6 + y]  += Axs[x]*Bxp[y];
	}

	int gy = iby+ty2;
	#pragma unroll
	for(int y=0; y<6; y++, gy+=16)
	{
		int gx = ibx+tx2;
		#pragma unroll
		for(int x=0; x<6; x++, gx+=16)
		{
			if (gx < m && gy < n)
				C[x*16] = alpha*Cb[y+x*6] + beta * C[x*16];
		}
		C+=ldc*16;
	}
}

//========================================================================

extern "C" __global__ void 
fermiSgemm_v2_kernel_TT(float *C, const float *A, const float *B, 
                        int m, int n,  int k,  int lda,  int ldb, 
                        int ldc, float alpha, float beta,
                        int offsetA, int offsetB) 
{
	const  int tx = threadIdx.x;
	const  int ty = threadIdx.y;

	const int iby = blockIdx.y * 96;
	const int ibx = blockIdx.x * 96;
	const int idt = ty * 64 + tx;

	const int tx2 = idt%16;
	const int ty2 = idt/16;

	__shared__ float Bb[16][97];
	__shared__ float Abs[96][17];

	float xxA[6];
	float xxB[6];

	int trackA = offsetA + __mul24( ibx + ty2, lda) + tx2;
	int trackB = offsetB + iby+ tx2 + __mul24(ty2, ldb);

	A += trackA; 
	B += trackB; 

	int tll = tx2; 

	#pragma unroll
	for(int y=0; y<6; y++)
		Abs[ty2+16*y][tx2] = /* (tll<k)* */ fetch_x_A(trackA +  lda*16*y, A);

	#pragma unroll
	for(int y=0; y<6; y++)
		Bb[ty2][tx2+16*y] = fetch_x_B(trackB+16*y, B);

	__syncthreads();

	const float *Bend = B + k*ldb - 16*ldb;

	float Axs[6];
	float Bxp[6];

	float Cb[36] = {0,0,0,0,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, 0,0,0,0,0,0,
		0,0,0,0,0,0, 0,0,0,0,0,0};
	do 
	{
		tll+=16;
		A += 16;
		B += 16*ldb;
		trackA+=16; 
		trackB+=16*ldb;
		
		#pragma unroll
		for( int y=0; y<6; y++)
			xxA[y] = /* (tll<k)* */ fetch_x_A(trackA + lda*y*16, A);

		#pragma unroll
		for( int y=0; y<6; y++)
			xxB[y] = fetch_x_B(trackB + 16*y, B);

		#pragma unroll 
		for( int j1=0;j1<16;j1++)
		{
			#pragma unroll
			for( int y=0; y<6; y++)
				Axs[y] =  Abs[tx2 + y*16][j1];

			#pragma unroll
			for( int y=0; y<6; y++)
				Bxp[y]= Bb[j1][ty2 + y*16];

			#pragma unroll 
			for( int x=0; x<6; x++)
				#pragma unroll 
				for( int y=0;y<6;y++)
					Cb[x*6+y] += Axs[x]*Bxp[y];
		}
		__syncthreads();

		#pragma unroll
		for( int y=0; y<6; y++)
			Abs[ty2 + 16*y][tx2] = xxA[y];

		#pragma unroll
		for( int y=0; y<6; y++)
			Bb[ty2][tx2+y*16] = xxB[y];

		__syncthreads();
	} 
	while (B < Bend);

	C += tx2 + ibx  + __mul24 (ty2 +  iby ,ldc);

	#pragma unroll 
	for( int j1=0; j1<16; j1++)
	{
		#pragma unroll
		for( int y=0; y<6; y++)
			Axs[y] =  Abs[tx2 + y*16][j1];

		#pragma unroll
		for( int y=0; y<6; y++)
			Bxp[y]= Bb[j1][ty2 + y*16];

		#pragma unroll 
		for( int x=0; x<6; x++)
			#pragma unroll 
			for( int y=0; y<6; y++)
				Cb[x*6+y]  += Axs[x]*Bxp[y];
	}

	int gy = iby + ty2;
	#pragma unroll
	for( int y=0; y<6; y++, gy+=16)
	{
		int gx = ibx + tx2; 
		#pragma unroll
		for(int x=0; x<6; x++, gx+=16)
		{
			if (gx < m && gy < n)
				C[x*16] = alpha*Cb[y+x*6] + beta * C[x*16];
		}

		C+=ldc*16;
	}
}


//========================================================================

extern "C" __global__ void 
fermiSgemm_v2_kernel_NT(float *C, const float *A, const float *B,  
                        int m, int n,  int k,  int lda,  int ldb,  
                        int ldc, float alpha, float beta,
                        int offsetA, int offsetB) 
{
	const  int tx = threadIdx.x;
	const  int ty = threadIdx.y;

	const int iby = blockIdx.y * 96;
	const int ibx = blockIdx.x * 96;
	const int idt = ty * 64 + tx;

	const int tx2= idt%16;
	const int ty2= idt/16;

	__shared__ float Bb[16][97];
	__shared__ float Abs[96][17];

	float xxA[6];
	float xxB[6];

	int trackA = offsetA + ibx +__mul24(ty2, lda) + tx2 ;
	int trackB = offsetB + iby + tx2 + __mul24(ty2, ldb);	
	
	A+= trackA; 
	B += trackB; 

	int tll = ty2;

	#pragma unroll
	for(int y=0; y<6; y++)
		Abs[tx2+ y*16][ty2] = /* (tll<k)* */ fetch_x_A(trackA + y*16, A);

	#pragma unroll
	for(int y=0; y<6; y++)
		Bb[ty2][tx2+16*y] = /* (tll<k)* */ fetch_x_B(trackB+16*y, B);

	__syncthreads();

	const float *Bend = B + k*ldb - 16*ldb;

	float Axs[6];
	float Bxp[6];

	float Cb[36] = {0,0,0,0,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, 0,0,0,0,0,0,
		0,0,0,0,0,0, 0,0,0,0,0,0};
	do 
	{
		tll += 16;
		A += lda *16  ;
		B += 16*ldb;
		trackA+=16*lda ; 
		trackB+=16*ldb;

		#pragma unroll
		for( int y=0; y<6; y++)
			xxA[y] = /* (tll<k)* */ fetch_x_A(trackA + y*16, A);	
			// tll same in the NN case

		#pragma unroll
		for( int y=0; y<6; y++)
			xxB[y] = /* (tll<k)* */ fetch_x_B( trackB + 16*y, B);

		#pragma unroll 
		for( int j1=0;j1<16;j1++)
		{
			#pragma unroll
			for( int y=0; y<6; y++)
				Bxp[y]= Bb[j1][ty2 + y*16];
			#pragma unroll
			for( int y=0; y<6; y++)
				Axs[y] =  Abs[tx2 + y*16][j1] ;

			#pragma unroll 
			for( int x=0; x<6; x++)
				#pragma unroll 
				for( int y=0; y<6; y++)
					Cb[x*6+y] += Axs[x]*Bxp[y];
		}
		__syncthreads();

		#pragma unroll
		for( int y=0; y<6; y++)
			Abs[tx2 + y*16][ty2] = xxA[y]; 

		#pragma unroll
		for( int y=0; y<6; y++)
			Bb[ty2][tx2+y*16] = xxB[y];

		__syncthreads();
	} 
	while (B < Bend);

	C += tx2 + ibx + __mul24(ty2 + iby ,ldc);

	#pragma unroll 
	for(int j1=0; j1<16; j1++)
	{
		#pragma unroll
		for( int y=0; y<6; y++)
			Bxp[y] = Bb[j1][ty2 + y*16];

		#pragma unroll
		for( int y=0; y<6; y++)
			Axs[y] =  Abs[tx2 + y*16][j1] ;

		#pragma unroll 
		for( int x=0; x<6; x++)
			#pragma unroll 
			for( int y=0;y<6;y++)
				Cb[x*6+y]  += Axs[x]*Bxp[y];
	}

	int gy = iby + ty2;
	#pragma unroll
	for( int y=0; y<6; y++, gy+=16)
	{
		int gx = ibx + tx2; 
		#pragma unroll
		for(int x=0; x<6; x++, gx+=16)
		{
			if (gx < m && gy < n)
				C[x*16] = alpha*Cb[y + x*6] + beta * C[x*16];
		}

		C+=ldc*16;
	}
}

//=================================================================================

extern "C" void
magmablas_sgemm_fermi( char TRANSA, char TRANSB, int m , int n , int k , 
                       float alpha, const float *A, int lda, 
                                    const float *B, int ldb, 
                       float beta,        float *C, int ldc ) 
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

   Purpose
   =======

   SGEMM  performs one of the matrix-matrix operations

      C := alpha*op( A )*op( B ) + beta*C,

   where  op( X ) is one of

      op( X ) = X   or   op( X ) = X',

   alpha and beta are scalars, and A, B and C are matrices, with op( A )
   an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

   Parameters
   ==========
   TRANSA - CHARACTER*1.
            On entry, TRANSA specifies the form of op( A ) to be used in
            the matrix multiplication as follows:
               TRANSA = 'N' or 'n',  op( A ) = A.
               TRANSA = 'T' or 't',  op( A ) = A'.
               TRANSA = 'C' or 'c',  op( A ) = A'.
            Unchanged on exit.

   TRANSB - CHARACTER*1.
            On entry, TRANSB specifies the form of op( B ) to be used in
            the matrix multiplication as follows:
               TRANSB = 'N' or 'n',  op( B ) = B.
               TRANSB = 'T' or 't',  op( B ) = B'.
               TRANSB = 'C' or 'c',  op( B ) = B'.
            Unchanged on exit.

   M      - INTEGER.
            On entry,  M  specifies  the number  of rows  of the  matrix
            op( A )  and of the  matrix  C.  M  must  be at least  zero.
            Unchanged on exit.

   N      - INTEGER.
            On entry,  N  specifies the number  of columns of the matrix
            op( B ) and the number of columns of the matrix C. N must be
            at least zero.
            Unchanged on exit.

   K      - INTEGER.
            On entry,  K  specifies  the number of columns of the matrix
            op( A ) and the number of rows of the matrix op( B ). K must
            be at least  zero.
            Unchanged on exit.

   ALPHA  - SINGLE PRECISION.
            On entry, ALPHA specifies the scalar alpha.
            Unchanged on exit.

   A      - SINGLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
            k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
            Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
            part of the array  A  must contain the matrix  A,  otherwise
            the leading  k by m  part of the array  A  must contain  the
            matrix A.
            Unchanged on exit.

   LDA    - INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program. When  TRANSA = 'N' or 'n' then
            LDA must be at least  max( 1, m ), otherwise  LDA must be at
            least  max( 1, k ).
            Unchanged on exit.

   B      - SINGLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
            n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
            Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
            part of the array  B  must contain the matrix  B,  otherwise
            the leading  n by k  part of the array  B  must contain  the
            matrix B.
            Unchanged on exit.

   LDB    - INTEGER.
            On entry, LDB specifies the first dimension of B as declared
            in the calling (sub) program. When  TRANSB = 'N' or 'n' then
            LDB must be at least  max( 1, k ), otherwise  LDB must be at
            least  max( 1, n ).
            Unchanged on exit.

   BETA   - SINGLE PRECISION.
            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
            supplied as zero then C need not be set on input.
            Unchanged on exit.

   C      - SINGLE PRECISION array of DIMENSION ( LDC, n ).
            Before entry, the leading  m by n  part of the array  C must
            contain the matrix  C,  except when  beta  is zero, in which
            case C need not be set on entry.
            On exit, the array  C  is overwritten by the  m by n  matrix
            ( alpha*op( A )*op( B ) + beta*C ).

   LDC    - INTEGER.
            On entry, LDC specifies the first dimension of C as declared
            in  the  calling  (sub)  program.   LDC  must  be  at  least
            max( 1, m ).
            Unchanged on exit.
   =====================================================================    */

	if (m<=0 || n<=0 || k<=0)
	   return;

	size_t offsetA = 0;
	size_t offsetB = 0;

	int TransA = 1, TransB = 1;
	if (TRANSA == 'N' ||  TRANSA == 'n')
	   TransA = 0;
	if (TRANSB == 'N' ||  TRANSB == 'n')
	   TransB = 0;

	size_t sizeA = (size_t) lda * (size_t) (!TransA ? k : m);
	size_t sizeB = (size_t) ldb * (size_t) (!TransB ? n : k);

 	size_t CUBLAS_MAX_1DBUF_SIZE = (1 << 27) - 512;
	if (sizeA>=CUBLAS_MAX_1DBUF_SIZE ||
			sizeB>=CUBLAS_MAX_1DBUF_SIZE )
	{
		//printf("Exceeding texuture limit (CUBLAS_MAX_1DBUF_SIZE=%ld), using cublasSgemm\n", CUBLAS_MAX_1DBUF_SIZE);
		cublasSgemm(TRANSA, TRANSB, m, n, k, alpha,
				A, lda, B, ldb,
				beta, C, ldc);
		return;
	}

	cudaError_t errt;
	errt = cudaBindTexture(&offsetA, tex_x_float_A, (int2 *)A, 
			sizeA * sizeof(A[0]));
	if( errt != cudaSuccess) printf("can not bind to texture \n");

	errt = cudaBindTexture(&offsetB, tex_x_float_B, (int2 *)B, 
			sizeB * sizeof(B[0]));
	if( errt != cudaSuccess) printf("can not bind to texture \n");

	dim3 threads( 64, 4 );
	dim3 grid(m/(96)+(m%(96)!=0),n/(96)+(n%(96)!=0));

	offsetA = offsetA/sizeof(A[0]);
	offsetB = offsetB/sizeof(B[0]);

	if ( TransB ) 
	   if( !TransA ) 
		fermiSgemm_v2_kernel_NT<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
					                    ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);
	   else
		fermiSgemm_v2_kernel_TT<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
					                    ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);
	else
	   if( !TransA ) 
		fermiSgemm_v2_kernel_NN<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
					                    ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);
	   else
		fermiSgemm_v2_kernel_TN<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
					                    ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);

	cudaUnbindTexture ( tex_x_float_A ) ;
	cudaUnbindTexture ( tex_x_float_B ) ;
}
