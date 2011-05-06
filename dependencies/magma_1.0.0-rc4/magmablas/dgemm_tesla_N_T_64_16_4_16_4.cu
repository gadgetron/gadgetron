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
dgemm_kernel_N_T_64_16_4_16_4(double *C, const double *A, const double *B, 
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
       C = alpha* A*B^T  + beta * C

    B is put into shared memory
    Parameters Used:
        blk_M=64 blk_N=16 blk_K=4 nthd_x=16 nthd_y=4

    This code should run for any matrix size.
    ===============================================================  */

	const int tx = threadIdx.x;
	const int ty = threadIdx.y;

	const int ibx = blockIdx.x * 64;
	const int iby = blockIdx.y *16;
	

	const int idt = ty * 16 + tx;

	if( iby + tx >=n )
		B+= iby+0;
	else
		B+= iby+tx;
	/*
		Taking care of boundary cases where K<4.
	*/
	if( ty >=k ) 
		B+= __mul24( 0,ldb);
	else
		B+= __mul24( ty,ldb);
		

	if( ibx + idt >= m ) 
		A += ibx + 0 ; 
	else
		A += ibx + idt;


	int s2=lda, s3=2*lda, s4=3*lda ; 

	switch (k){
		case 1: 
			s2=0; s3=0;s4=0 ; 
			break ; 
		case 2:
			s2=lda; s3=0;s4=0 ; 
			break ; 
		case 3:  
			s2=lda; s3=2*lda;s4=0 ; 
			break ;
	}
	
	C += ibx +idt +__mul24( iby,ldc);

	double Ap[4]={A[0], A[s2], A[s3], A[s4]};

	double b=B[0];

	const double *Bend = B + ldb*(k-k%4);

	B+=4*ldb;
	A+=4*lda;

	__shared__ double Bb[4][16];

	double Cb[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	if(k>7)
	do {

		double Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};

		Bb[ty][tx]=b;

		__syncthreads();

		Ap[0] = A[0];
		Ap[1] = A[s2];
		Ap[2] = A[s3];
		Ap[3] = A[s4];

		 b=B[0];


		daxpy(Ab[0], &Bb[0][0], Cb);
		daxpy(Ab[1], &Bb[1][0], Cb);
		daxpy(Ab[2], &Bb[2][0], Cb);
		daxpy(Ab[3], &Bb[3][0], Cb);

		A+=4*lda;
		B += 4*ldb;

		__syncthreads();
	} while (B < Bend);

	if(k>3){

		Bb[ty][tx]=b;
		int k1 = k-k%4;

		if( (k1+ty) >=k)
		    B-=4*ldb;
		else 
		    B-=0*ldb;

		if( (k1+0) >= k ) {s2=0  ;s3=0*lda;s4=0;A-=4*lda;} else
		if( (k1+1) >= k ) {s2=0  ;s3=0*lda;s4=0;A-=0*lda;} else
		if( (k1+2) >= k ) {s2=lda;s3=0*lda;s4=0;A-=0*lda;} else
		if( (k1+3) >= k ) {s2=lda;s3=2*lda;s4=0;A-=0*lda;} 
			
		__syncthreads();


	         b=B[0];

		daxpy(Ap[0], &Bb[0][0], Cb);        Ap[0] = A[0];
		daxpy(Ap[1], &Bb[1][0], Cb);        Ap[1] = A[s2];
		daxpy(Ap[2], &Bb[2][0], Cb);        Ap[2] = A[s3];
		daxpy(Ap[3], &Bb[3][0], Cb);        Ap[3] = A[s4];
	
	}

	k=k%4;

	if ( k!=0){

		__syncthreads();

		Bb[ty][tx]=b;

		__syncthreads();

		for(int i=0;i<k;i++){
			daxpy(Ap[i],&Bb[i][0], Cb);
		}
	}



	if( (iby+16)>=n) { 
		lda = n-iby;
	}
	else{
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
magmablas_dgemm_kernel_N_T_64_16_4_16_4(double *C,
                                        const double *A,
                                        const double *B,
                                        int m, int n, int k, 
                                        int lda, int ldb, int ldc,
                                        double alpha, double beta)
{	
        dim3 threads( 16, 4 );
        dim3 grid(m/64+(m%64!=0),n/16+(n%16!=0));
        dgemm_kernel_N_T_64_16_4_16_4<<< grid, threads >>>(C, A, B, 
                                                           m, n, k, 
                                                           lda, ldb, ldc, 
                                                           alpha, beta);
}

