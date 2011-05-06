/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/

/*
    blk_M=64 blk_N=64 blk_K=16 nthd_x=64 nthd_y=4
*/
#include "common_magma.h"
#include "commonblas_d.h"

#define magmablas_dgemm_fermi magmablas_dgemm

texture<int2,1>  tex_x_double_A;
texture<int2,1>  tex_x_double_B;

static __inline__ __device__ double fetch_x_A(const int& i)
{
  register int2  v = tex1Dfetch(tex_x_double_A, i);
  return __hiloint2double(v.y, v.x);
}

static __inline__ __device__ double fetch_x_B(const int& i)
{
  register int2  v = tex1Dfetch(tex_x_double_B, i);
  return __hiloint2double(v.y, v.x);
}

extern "C" __global__ void 
fermiDgemm_v2_kernel_NN(double *C, const double *A, const double *B,  
                        int m, int n, int k, int lda, int ldb,  
                        int ldc, double alpha, double beta,
                        int offsetA, int offsetB) 
{
	const  int tx = threadIdx.x;
	const  int ty = threadIdx.y;

	const int iby = blockIdx.y * 64;
	const int ibx = blockIdx.x * 64;
	const int idt = ty * 64 + tx;

	const int tx2 = idt%16;
	const int ty2 = idt/16;

	__shared__ double Abs[64][17];
	__shared__ double  Bb[16][65];

	int tll = ty2;
	double xxA[4];
	double xxB[4];

	int trackA = offsetA + ibx +__mul24( ty2, lda) + tx2 ;
	A += trackA; 

	int trackB = offsetB + tx2+ __mul24(iby + ty2 * 4, ldb );
	B += trackB;

	#pragma unroll
	for(int y=0; y<4; y++)
		Abs[tx2+ y*16][ty2] = /* (tll<k)* */ fetch_x_A(trackA + y*16) ;

	#pragma unroll
	for(int y=0; y<4; y++)
		Bb[tx2][ty2*4+y] = fetch_x_B( trackB + y * ldb) ;

	__syncthreads();

	double Axs[4];
	double Bxp[4];

	double Cb[16] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

	int k1;
	for(k1=0; k1<(k-16); k1+=16)
	{
		tll+=16;
		A += lda *16  ;
		B += 16;
		trackA += 16*lda ; 
		trackB += 16;

		#pragma unroll
		for( int y=0; y<4; y++)
			xxA[y] = /* (tll<k)* */ fetch_x_A(trackA + y*16);

		#pragma unroll
		for( int y=0; y<4; y++)
			xxB[y] = fetch_x_B( trackB + y*ldb);

		#pragma unroll 
		for( int j1=0;j1<16;j1++)
		{
			#pragma unroll
			for( int y=0; y<4; y++)
				Axs[y] =  Abs[tx2+y*16][j1] ;

			#pragma unroll
			for( int y=0; y<4; y++)
				Bxp[y]= Bb[j1][ty2+y*16];


			#pragma unroll 
			for( int x=0; x<4; x++)
			{
				#pragma unroll 
				for( int y=0; y<4; y++)
				{
					Cb[x*4+y]  += Axs[x]*Bxp[y];
				}
			}
		}

		__syncthreads();
		
		#pragma unroll
		for(int y=0; y<4; y++)
			Abs[tx2+y*16][ty2] = xxA[y]; 

		#pragma unroll
		for(int y=0; y<4; y++)
			Bb[tx2][ty2*4 + y] = xxB[y];

		__syncthreads();
	}

	C += tx2 + ibx  + __mul24 (ty2 +  iby ,ldc);

	#pragma unroll 
	for(int j1=0;j1<16;j1++)
	{

		#pragma unroll
		for( int y=0; y<4; y++)
			Axs[y] =  Abs[tx2 + y*16][j1] ;

		#pragma unroll
		for( int y=0; y<4; y++)
			Bxp[y]= Bb[j1][ty2 + y*16];

		#pragma unroll 
		for( int x=0; x<4; x++)
		{
			#pragma unroll 
			for( int y=0;y<4; y++)
			{
				Cb[x*4 + y]  += Axs[x]*Bxp[y];
			}
		}
	}

	int gy = iby + ty2;
	#pragma unroll
	for( int y=0;y<4;y++, gy+=16)
	{
		int gx = ibx + tx2; 
		#pragma unroll
		for(int x=0;x<4;x++, gx+=16)
		{
			if (gx < m && gy < n)
				C[x*16] = alpha*Cb[y+x*4] + beta * C[x*16];
		}

		C += ldc*16;
	}
}

extern "C" __global__ void 
fermiDgemm_v2_kernel_TN(double *C, const double *A, const double *B,
                        int m, int n,  int k,  int lda,  int ldb,  
                        int ldc, double alpha, double beta,
                        int offsetA, int offsetB) 
{
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;

    const int iby = blockIdx.y * 64;
    const int ibx = blockIdx.x * 64;
    const int idt = ty * 64 + tx;

    const int tx2 = idt%16;
    const int ty2 = idt/16;

    __shared__ double Bb[16][65];
    __shared__ double Abs[64][17];

    double xxA[4];
    double xxB[4];

    int trackA = offsetA + tx2 + __mul24( ibx + ty2*4, lda ); 
    int trackB = offsetB + tx2 + __mul24( iby + ty2*4, ldb ); 

    A+= trackA; 
    B+= trackB; 

	int tll = tx2;

    #pragma unroll
    for(int y=0; y<4; y++)
		Abs[ty2*4+y][tx2] =  (tll<k)*  fetch_x_A(trackA + y*lda);

    #pragma unroll
    for(int y=0; y<4; y++)
		Bb[tx2][ty2*4+y] = /* (tll<k)* */ fetch_x_B( trackB + y*ldb );

    __syncthreads();
   
    double Axs[4];
    double Bxp[4];

    double Cb[16] = {0,0,0,0,    0,0,0,0, 0,0,0,0, 0,0,0,0};

    int k1;
    for(k1=0; k1<(k-16); k1+=16)
    {
		tll +=16;
		B += 16;
		A += 16  ;
		trackA+=16 ; 
		trackB+=16;

		#pragma unroll
		for(int y=0; y<4; y++)
			xxA[y] =  (tll<k)*  fetch_x_A(trackA + y*lda);

		#pragma unroll
		for(int y=0; y<4; y++)
			xxB[y] = /* (tll<k)* */ fetch_x_B(trackB + y*ldb);

		#pragma unroll 
		for(int j1=0;j1<16;j1++)
		{
			#pragma unroll
			for(int y=0; y<4; y++)
				Axs[y] =  Abs[tx2+y*16][j1];

			#pragma unroll
			for(int y=0; y<4; y++)
				Bxp[y]= Bb[j1][ty2+y*16];

			#pragma unroll 
			for(int x=0; x<4; x++)
			{
				#pragma unroll 
				for(int y=0; y<4; y++)
				{
					Cb[x*4+y]  += Axs[x]*Bxp[y];
				}
			}
		}
		__syncthreads();

		#pragma unroll
		for(int y=0; y<4; y++)
			Abs[ty2*4+y][tx2] = xxA[y];

		#pragma unroll
		for(int y=0; y<4; y++)
			Bb[tx2][ty2*4+y] =xxB[y];
		__syncthreads();
	}

	C += tx2 + ibx  + __mul24 (ty2 + iby ,ldc);

	#pragma unroll 
	for(int j1=0; j1<16; j1++)
	{
		#pragma unroll
		for(int y=0; y<4; y++)
			Axs[y] = Abs[tx2+y*16][j1];

		#pragma unroll
		for(int y=0; y<4; y++)
			Bxp[y]= Bb[j1][ty2+y*16];

		#pragma unroll 
		for(int x=0; x<4; x++)
		{
			#pragma unroll 
			for(int y=0; y<4; y++)
			{
				Cb[x*4+y] += Axs[x]*Bxp[y];
			}
		}
	}

	int gy = iby+ty2;
	#pragma unroll
	for(int y=0;y<4;y++, gy+=16)
	{
		int gx = ibx+tx2;
		#pragma unroll
		for(int x=0;x<4;x++, gx+=16)
		{
			if (gx < m && gy < n)
			C[x*16] =alpha*Cb[y+x*4] + beta * C[x*16];
		}
		C+=ldc*16;
	}
}

extern "C" __global__ void 
fermiDgemm_v2_kernel_TT(double *C, const double *A, const double *B, 
                        int m, int n,  int k,  int lda,  int ldb, 
                        int ldc, double alpha, double beta,
                        int offsetA, int offsetB) 
{
	const  int tx = threadIdx.x;
	const  int ty = threadIdx.y;

	const int iby = blockIdx.y * 64;
	const int ibx = blockIdx.x * 64;
	const int idt = ty * 64 + tx;

	const int tx2 = idt%16;
	const int ty2 = idt/16;

	__shared__ double Bb[16][65];
	__shared__ double Abs[64][17];

	double xxA[4];
	double xxB[4];

	int trackA = offsetA + __mul24( ibx + ty2, lda) + tx2;
	int trackB = offsetB + iby+ tx2 + __mul24(ty2, ldb);

	A += trackA; 
	B += trackB; 

	int tll = tx2; 

	#pragma unroll
	for(int y=0; y<4; y++)
		Abs[ty2+16*y][tx2] = /* (tll<k)* */ fetch_x_A(trackA +  lda*16*y);

	#pragma unroll
	for(int y=0; y<4; y++)
		Bb[ty2][tx2+16*y] = fetch_x_B(trackB+16*y);

	__syncthreads();

	double Axs[4];
	double Bxp[4];

	double Cb[16] = {0,0,0,0, 0,0,0,0,  0,0,0,0, 0,0,0,0};
	
        int k1;
	for(k1=0; k1<(k-16); k1+=16)
	{
		tll+=16;
		A += 16;
		B += 16*ldb;
		trackA+=16; 
		trackB+=16*ldb;
		
		#pragma unroll
		for( int y=0; y<4; y++)
			xxA[y] = /* (tll<k)* */ fetch_x_A(trackA + lda*y*16);

		#pragma unroll
		for( int y=0; y<4; y++)
			xxB[y] = fetch_x_B(trackB + 16*y);

		#pragma unroll 
		for( int j1=0;j1<16;j1++)
		{
			#pragma unroll
			for( int y=0; y<4; y++)
				Axs[y] =  Abs[tx2 + y*16][j1];

			#pragma unroll
			for( int y=0; y<4; y++)
				Bxp[y]= Bb[j1][ty2 + y*16];

			#pragma unroll 
			for( int x=0; x<4; x++)
				#pragma unroll 
				for( int y=0;y<4;y++)
					Cb[x*4+y] += Axs[x]*Bxp[y];
		}
		__syncthreads();

		#pragma unroll
		for( int y=0; y<4; y++)
			Abs[ty2 + 16*y][tx2] = xxA[y];

		#pragma unroll
		for( int y=0; y<4; y++)
			Bb[ty2][tx2+y*16] = xxB[y];

		__syncthreads();
	} 

	C += tx2 + ibx  + __mul24 (ty2 +  iby ,ldc);

	#pragma unroll 
	for( int j1=0; j1<16; j1++)
	{
		#pragma unroll
		for( int y=0; y<4; y++)
			Axs[y] =  Abs[tx2 + y*16][j1];

		#pragma unroll
		for( int y=0; y<4; y++)
			Bxp[y]= Bb[j1][ty2 + y*16];

		#pragma unroll 
		for( int x=0; x<4; x++)
			#pragma unroll 
			for( int y=0; y<4; y++)
				Cb[x*4+y]  += Axs[x]*Bxp[y];
	}

	int gy = iby + ty2;
	#pragma unroll
	for( int y=0; y<4; y++, gy+=16)
	{
		int gx = ibx + tx2; 
		#pragma unroll
		for(int x=0; x<4; x++, gx+=16)
		{
			if (gx < m && gy < n)
				C[x*16] = alpha*Cb[y+x*4] + beta * C[x*16];
		}

		C+=ldc*16;
	}
}

	
extern "C" __global__ void 
fermiDgemm_v2_kernel_NT(double *C, const double *A, const double *B,  
                        int m, int n,  int k,  int lda,  int ldb,  
                        int ldc, double alpha, double beta,
                        int offsetA, int offsetB) 
{
	const  int tx = threadIdx.x;
	const  int ty = threadIdx.y;

	const int iby = blockIdx.y * 64;
	const int ibx = blockIdx.x * 64;
	const int idt = ty * 64 + tx;

	const int tx2= idt%16;
	const int ty2= idt/16;

	__shared__ double Bb[16][65];
	__shared__ double Abs[64][17];

	double xxA[4];
	double xxB[4];

	int trackA = offsetA + ibx +__mul24(ty2, lda) + tx2 ;
	int trackB = offsetB + iby + tx2 + __mul24(ty2, ldb);	
	
	A+= trackA; 
	B += trackB; 

	int tll = ty2;

	#pragma unroll
	for(int y=0; y<4; y++)
		Abs[tx2+ y*16][ty2] = /* (tll<k)* */ fetch_x_A(trackA + y*16);

	#pragma unroll
	for(int y=0; y<4; y++)
		Bb[ty2][tx2+16*y] = /* (tll<k)* */ fetch_x_B(trackB+16*y);

	__syncthreads();

	double Axs[4];
	double Bxp[4];

	double Cb[16] = {0,0,0,0, 0,0,0,0,  0,0,0,0, 0,0,0,0};
	
        int k1;
	for(k1=0; k1<(k-16); k1+=16)
	{
		tll += 16;
		A += lda *16  ;
		B += 16*ldb;
		trackA+=16*lda ; 
		trackB+=16*ldb;

		#pragma unroll
		for( int y=0; y<4; y++)
			xxA[y] = /* (tll<k)* */ fetch_x_A(trackA + y*16);

		#pragma unroll
		for( int y=0; y<4; y++)
			xxB[y] = /* (tll<k)* */ fetch_x_B( trackB + 16*y);

		#pragma unroll 
		for( int j1=0;j1<16;j1++)
		{
			#pragma unroll
			for( int y=0; y<4; y++)
				Bxp[y]= Bb[j1][ty2 + y*16];
			#pragma unroll
			for( int y=0; y<4; y++)
				Axs[y] =  Abs[tx2 + y*16][j1] ;

			#pragma unroll 
			for( int x=0; x<4; x++)
				#pragma unroll 
				for( int y=0; y<4; y++)
					Cb[x*4+y] += Axs[x]*Bxp[y];
		}
		__syncthreads();

		#pragma unroll
		for( int y=0; y<4; y++)
			Abs[tx2 + y*16][ty2] = xxA[y]; 

		#pragma unroll
		for( int y=0; y<4; y++)
			Bb[ty2][tx2+y*16] = xxB[y];

		__syncthreads();
	} 

	C += tx2 + ibx + __mul24(ty2 + iby ,ldc);

	#pragma unroll 
	for(int j1=0; j1<16; j1++)
	{
		#pragma unroll
		for( int y=0; y<4; y++)
			Bxp[y] = Bb[j1][ty2 + y*16];

		#pragma unroll
		for( int y=0; y<4; y++)
			Axs[y] =  Abs[tx2 + y*16][j1] ;

		#pragma unroll 
		for( int x=0; x<4; x++)
			#pragma unroll 
			for( int y=0;y<4;y++)
				Cb[x*4+y]  += Axs[x]*Bxp[y];
	}

	int gy = iby + ty2;
	#pragma unroll
	for( int y=0; y<4; y++, gy+=16)
	{
		int gx = ibx + tx2; 
		#pragma unroll
		for(int x=0; x<4; x++, gx+=16)
		{
			if (gx < m && gy < n)
				C[x*16] = alpha*Cb[y + x*4] + beta * C[x*16];
		}

		C+=ldc*16;
	}
}

extern "C" void
magmablas_dgemm_fermi( char TRANSA, char TRANSB, int m , int n , int k , 
                       double alpha, const double *A, int lda, 
                                     const double *B, int ldb, 
                       double beta,        double *C, int ldc ) 
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

   Purpose
   =======

   DGEMM  performs one of the matrix-matrix operations

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

   ALPHA  - DOUBLE PRECISION.
            On entry, ALPHA specifies the scalar alpha.
            Unchanged on exit.

   A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
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

   B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
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

   BETA   - DOUBLE PRECISION.
            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
            supplied as zero then C need not be set on input.
            Unchanged on exit.

   C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
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
 	
        // size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512) / 2;
        size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512);
	if (sizeA>=CUBLAS_MAX_1DBUF_SIZE ||
			sizeB>=CUBLAS_MAX_1DBUF_SIZE )
	{
	//	printf("Exceeding texuture limit (CUBLAS_MAX_1DBUF_SIZE=%ld), using cublasSgemm\n", CUBLAS_MAX_1DBUF_SIZE);
		cublasDgemm(TRANSA, TRANSB, m, n, k, alpha,
				A, lda, B, ldb,
				beta, C, ldc);
		return;
	}

	cudaError_t  errt;
	errt = cudaBindTexture(&offsetA, tex_x_double_A, (int2 *)A,
                               sizeA * sizeof(A[0]));
	if( errt != cudaSuccess) 
	{
		printf("can not bind to texture \n");
		return;
	}

	errt = cudaBindTexture(&offsetB, tex_x_double_B, (int2 *)B,
                               sizeB * sizeof(B[0]));
	if( errt != cudaSuccess)
	{
		printf("can not bind to texture \n");
		return;
	}

	dim3 threads( 64, 4 );
	dim3 grid(m/(64)+(m%(64)!=0),n/(64)+(n%(64)!=0));

	offsetA = offsetA/sizeof(A[0]);
        offsetB = offsetB/sizeof(B[0]);

	if ( TransB ) 
	   if ( !TransA ) 
		fermiDgemm_v2_kernel_NT<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
					                    ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);
	   else
		fermiDgemm_v2_kernel_TT<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
		                                            ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);
	else
	   if ( !TransA ) 
		fermiDgemm_v2_kernel_NN<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
                                                            ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);
	   else
		fermiDgemm_v2_kernel_TN<<< grid, threads>>>(C, A, B, m, n, k, lda, ldb, 
                                                            ldc, alpha, beta,
                                                            (int)offsetA, (int)offsetB);

	cudaUnbindTexture ( tex_x_double_A ) ;
	cudaUnbindTexture ( tex_x_double_B ) ;
}
