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
dgemm_kernel_T_N_32_32_8_8_8(double *C, const double *A, const double *B,
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
       C = alpha* A^T*B  + beta * C

    B is put into shared memory
    Parameters Used:
        blk_M=32 blk_N=32 blk_K=8 nthd_x=8 nthd_y=8

    This code should run for any matrix size.
    ===============================================================  */

	const int ibx = blockIdx.x *32;
	const int iby = blockIdx.y *32;

	const int tx =  threadIdx.y ;
	const int ty =  threadIdx.x ;


	int idt = tx * 8 + ty;

	if( ty >=k ) 		
		A += __mul24(ibx ,lda)+0;
	else
		A += __mul24(ibx ,lda)+ty;
	
	if( (ibx + tx ) >= m ) 
		A += __mul24(0,lda);
	else 
		A += __mul24(tx,lda);

	if( (iby+tx) >=n ) 
		B+= __mul24(iby+0,ldb);
	else
		B+= __mul24(iby+tx,ldb) ;
	if( ty>=k)
		B+=0;
	else
		B+= ty;

	C += ibx +idt%32 +__mul24( iby+16*(idt/32),ldc);

	lda = lda *8 ;
	ldb = ldb *8 ;


	int as1=0, as2=lda, as3=2*lda , as4 =3*lda; 
	int bs1=0 , bs2=ldb , bs3=2*ldb , bs4=3*ldb ;


	switch(k){
		case 1: as2=0   ; as3 = 0*lda;as4=0; bs2=0   ; bs3 = 0*ldb; bs4=0; break;
		case 2: as2=lda ; as3 = 0*lda;as4=0; bs2=ldb ; bs3 = 0*ldb; bs4=0; break;
		case 3: as2=lda ; as3 = 2*lda;as4=0; bs2=ldb ; bs3 = 2*ldb; bs4=0; break;
	}

	if( (ibx + tx    ) >=m ) { as1=0;   as2=0*lda; as3=0*lda ; as4 =0*lda; } else
	if( (ibx + tx +8 ) >=m ) { as1=0;   as2=0*lda; as3=0*lda ; as4 =0*lda; } else
	if( (ibx + tx +16) >=m ) { as1=0;   as2=1*lda; as3=0*lda ; as4 =0*lda; } else
	if( (ibx + tx +24) >=m ) { as1=0;   as2=1*lda; as3=2*lda ; as4 =0*lda; } 
		

	if( (iby + tx    ) >=n ) { bs1=0;   bs2=0*ldb; bs3=0*ldb ; bs4 =0*ldb; } else
	if( (iby + tx +8 ) >=n ) { bs1=0;   bs2=0*ldb; bs3=0*ldb ; bs4 =0*ldb; } else
	if( (iby + tx +16) >=n ) { bs1=0;   bs2=1*ldb; bs3=0*ldb ; bs4 =0*ldb; } else
	if( (iby + tx +24) >=n ) { bs1=0;   bs2=1*ldb; bs3=2*ldb ; bs4 =0*ldb; } 

 
	double b= B[bs1];
	double b1=B[bs2];
	double b2=B[bs3];
	double b3=B[bs4];


	double Ap[4]={A[as1], A[as2], A[as3],A[as4]};

	const double *Bend = B + (k-k%8);

	B+=8;
	A+=8;

	__shared__ double Bb[8][33];
	__shared__ double ABb[32][9];
	
	double Cb[16] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
        	

	const int l = 17*(idt/32) ;
	int idt1 = idt ;	
	idt = idt % 32 ; 
	if(k>15){
	do {
		

		Bb[ty][tx   ] = b;
		Bb[ty][tx+8 ] = b1;
		Bb[ty][tx+17] = b2;
		Bb[ty][tx+25] = b3;

		ABb[tx   ][ty] = Ap[0];
		ABb[tx+8 ][ty] = Ap[1];
		ABb[tx+16][ty] = Ap[2];
		ABb[tx+24][ty] = Ap[3];
	

		__syncthreads();
		
		
		daxpy(ABb[idt][0], &Bb[0][l], Cb);Ap[0]=A[as1]; 
		daxpy(ABb[idt][1], &Bb[1][l], Cb);Ap[1]=A[as2];
		daxpy(ABb[idt][2], &Bb[2][l], Cb);Ap[2]=A[as3];
		daxpy(ABb[idt][3], &Bb[3][l], Cb);Ap[3]=A[as4]; 
	
		daxpy(ABb[idt][4], &Bb[4][l], Cb);
		b=B[bs1];
		daxpy(ABb[idt][5], &Bb[5][l], Cb);
		b1=B[bs2];
		daxpy(ABb[idt][6], &Bb[6][l], Cb); 
		b2=B[bs3];
		daxpy(ABb[idt][7], &Bb[7][l], Cb); 
		b3=B[bs4];

		B += 8;
		A += 8;

		__syncthreads();

	} while (B < Bend);
	}
	 if(k>7){

		Bb[ty][tx   ] = b;
		Bb[ty][tx+8 ] = b1;
		Bb[ty][tx+17] = b2;
		Bb[ty][tx+25] = b3;

		ABb[tx   ][ty] = Ap[0];
		ABb[tx+8 ][ty] = Ap[1];
		ABb[tx+16][ty] = Ap[2];
		ABb[tx+24][ty] = Ap[3];

		__syncthreads();
		as1 = k-k%8;
		
		if(as1+ty>=k){ bs1=0*ldb;bs2=0*ldb;bs3=0*ldb;bs4=0*ldb;B-=8;}		
		if(as1+ty>=k){ as1=0*lda;as2=0*lda;as3=0*lda;as4=0*lda;A-=8;}		

		as1=0;
		daxpy(ABb[idt][0], &Bb[0][l], Cb);
		Ap[0]=A[as1]; 
		daxpy(ABb[idt][1], &Bb[1][l], Cb);
		Ap[1]=A[as2];
		daxpy(ABb[idt][2], &Bb[2][l], Cb);
		Ap[2]=A[as3];
		daxpy(ABb[idt][3], &Bb[3][l], Cb);
		Ap[3]=A[as4]; 
	
		daxpy(ABb[idt][4], &Bb[4][l], Cb);
		b=B[bs1];
		daxpy(ABb[idt][5], &Bb[5][l], Cb);
		b1=B[bs2];
		daxpy(ABb[idt][6], &Bb[6][l], Cb); 
		b2=B[bs3];
		daxpy(ABb[idt][7], &Bb[7][l], Cb); 
		b3=B[bs4];
	}
	k=k%8;
	if ( k!=0){
		__syncthreads();

		Bb[ty][tx]= b;
		Bb[ty][tx+8] = b1;
		Bb[ty][tx+17] = b2;
		Bb[ty][tx+25] = b3;

		ABb[tx][ty]= Ap[0];
		ABb[tx+8][ty] = Ap[1];
		ABb[tx+16][ty] = Ap[2];
		ABb[tx+24][ty] = Ap[3];
		__syncthreads();

		for(int i=0;i<k;i++){
			daxpy(ABb[idt][i],&Bb[i][l], Cb);
		}
	}

	if( (iby+16*(idt1/32+1))>=n) { 
		lda = n-iby-16*(idt1/32);
	}
	else    {
		lda = 16;
	}
	if( (ibx+idt) >= m )
		lda = 0 ;
	else lda = lda ;
	
	switch(lda){
		case 16:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			C[9*ldc] =alpha*Cb[9] + beta * C[9*ldc];
			C[10*ldc] =alpha*Cb[10] + beta * C[10*ldc];
			C[11*ldc] =alpha*Cb[11] + beta * C[11*ldc];
			C[12*ldc] =alpha*Cb[12] + beta * C[12*ldc];
			C[13*ldc] =alpha*Cb[13] + beta * C[13*ldc];
			C[14*ldc] =alpha*Cb[14] + beta * C[14*ldc];
			C[15*ldc] =alpha*Cb[15] + beta * C[15*ldc];
			break;
		case 15:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			C[9*ldc] =alpha*Cb[9] + beta * C[9*ldc];
			C[10*ldc] =alpha*Cb[10] + beta * C[10*ldc];
			C[11*ldc] =alpha*Cb[11] + beta * C[11*ldc];
			C[12*ldc] =alpha*Cb[12] + beta * C[12*ldc];
			C[13*ldc] =alpha*Cb[13] + beta * C[13*ldc];
			C[14*ldc] =alpha*Cb[14] + beta * C[14*ldc];
			break;
		case 14:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			C[9*ldc] =alpha*Cb[9] + beta * C[9*ldc];
			C[10*ldc] =alpha*Cb[10] + beta * C[10*ldc];
			C[11*ldc] =alpha*Cb[11] + beta * C[11*ldc];
			C[12*ldc] =alpha*Cb[12] + beta * C[12*ldc];
			C[13*ldc] =alpha*Cb[13] + beta * C[13*ldc];
			break;
		case 13:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			C[9*ldc] =alpha*Cb[9] + beta * C[9*ldc];
			C[10*ldc] =alpha*Cb[10] + beta * C[10*ldc];
			C[11*ldc] =alpha*Cb[11] + beta * C[11*ldc];
			C[12*ldc] =alpha*Cb[12] + beta * C[12*ldc];
			break; 
		case 12:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			C[9*ldc] =alpha*Cb[9] + beta * C[9*ldc];
			C[10*ldc] =alpha*Cb[10] + beta * C[10*ldc];
			C[11*ldc] =alpha*Cb[11] + beta * C[11*ldc];
			break;
		case 11:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			C[9*ldc] =alpha*Cb[9] + beta * C[9*ldc];
			C[10*ldc] =alpha*Cb[10] + beta * C[10*ldc];
			break;
		case 10:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			C[9*ldc] =alpha*Cb[9] + beta * C[9*ldc];
			break;
		case 9:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			C[8*ldc] =alpha*Cb[8] + beta * C[8*ldc];
			break;
		case 8:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			C[7*ldc] =alpha*Cb[7] + beta * C[7*ldc];
			break;
		case 7:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			C[6*ldc] =alpha*Cb[6] + beta * C[6*ldc];
			break;
		case 6:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			C[5*ldc] =alpha*Cb[5] + beta * C[5*ldc];
			break;
		case 5:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			C[4*ldc] =alpha*Cb[4] + beta * C[4*ldc];
			break;
		case 4:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			C[3*ldc] =alpha*Cb[3] + beta * C[3*ldc];
			break;
		case 3:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			C[2*ldc] =alpha*Cb[2] + beta * C[2*ldc];
			break;
		case 2:
			C[0] =alpha*Cb[0] + beta * C[0];
			C[1*ldc] =alpha*Cb[1] + beta * C[1*ldc];
			break;
		case 1:
			C[0] =alpha*Cb[0] + beta * C[0];
			break;
		case 0:
			break;
	}

}

extern "C" void
magmablas_dgemm_kernel_T_N_32_32_8_8_8(double *C, 
                                       const double *A,
                                       const double *B,
                                       int m, int n, int k,
                                       int lda, int ldb, int ldc, 
                                       double alpha, double beta)
{	
        dim3 threads( 8, 8 );
        dim3 grid(m/32+(m%32!=0),n/32+(n%32!=0));
        dgemm_kernel_T_N_32_32_8_8_8<<< grid, threads >>>(C, A, B, 
                                                          m, n, k, 
                                                          lda, ldb, ldc, 
                                                          alpha, beta);
}
