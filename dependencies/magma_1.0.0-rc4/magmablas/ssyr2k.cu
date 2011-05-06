/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"

#define block_M 32
#define block_N 32
#define thread_x 32
#define thread_y 2
#define unroll_f 16

__device__ void saxpy(float a,float *b, float *c) {
	#pragma unroll 
         for (int i = 0; i < unroll_f; i++) 
	      c[i] += a * b[i];
}

extern "C" __global__ void 
Ssyr2k_v16_ts_even_generic(float *C, const float *A, const float *B, 
                           int m, int in, int k, 
                           int lda, int ldb, int ldc, 
                           float alpha, float beta) 
{
	int tx = threadIdx.x;
	const int ty = threadIdx.y;
        int ibx = blockIdx.x ;
        int iby = blockIdx.y ;
        iby = (iby+ibx+3 ) % gridDim.y ;
        const int minor = iby&1  ;
        const bool  bottom = ibx >iby ;
        ibx = ( bottom ) ? (ibx-1) : ( iby + gridDim.y  );
        iby = ( bottom ) ?  iby    : ( blockIdx.x +minor + gridDim.y );
        if( iby > ibx ) iby = in ;
        ibx =ibx * block_M;
        iby =iby * block_N;

        const float *A1 = A ; 
        const float *B1 = B ; 
{ 
        B+= iby+tx;
        B+= __mul24( ty,ldb);
	A+= ibx + tx ; 
	C += ibx +tx +__mul24( iby+ty* unroll_f,ldc);

	float Ap[4];
	Ap[0]=A[0] ;
	Ap[1]=A[lda] ;
	Ap[2]=A[2*lda] ;
	Ap[3]=A[3*lda] ;

	float b=B[0];
        float b2=B[2*ldb];
	const float *Bend = B + ldb*k ;
	B+=4*ldb;
	A+=4*lda;
	__shared__ float Bb[4][block_N];
	float Cb[unroll_f] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		 __syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B+= 4*ldb;
		 __syncthreads();
	} while (B < Bend);
	 Bb[ty][tx]=b;
	 Bb[ty+2][tx]=b2;
	 __syncthreads();
	 saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);      
	 saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);      
	 saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);      
	 saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      
	 __syncthreads();
// -- 2nd Half            
        B=A1;
        A=B1;
        int tlda = lda ; lda = ldb ; ldb = tlda ; 
	B+= iby+tx;
	B+= __mul24( ty,ldb);
	A+= ibx + tx ; 
	Ap[0]=A[0] ;
	Ap[1]=A[lda] ;
	Ap[2]=A[2*lda] ;
	Ap[3]=A[3*lda] ;
	b=B[0];
        b2=B[2*ldb];
	const float *Bend1 = B + ldb*k;
	B+=4*ldb;
	A+=4*lda;
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		 __syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B += 4*ldb;
		 __syncthreads();
	} while (B < Bend1);
	Bb[ty][tx]=b;
	Bb[ty+2][tx]=b2;
	__syncthreads();
	saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);      
	saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);      
	saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);      
	saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      

	lda  = 0 ; 

	if ( iby < ibx ){
		tx = 15 ; 
	}
	else{
    		if ( tx > 15 ) {
        		if ( ty == 0 ){
		 		lda =  1 ; 
		  		tx=15;
			}
        	else{ 
           		lda = 1 ; 
	   		tx-=16;
        	}
    	}
    	else{
        	if ( ty  == 0 ) {
			lda = 1 ; 
		}
         	else {
			lda = 2 ;
			tx = 32  ;  
		}	
    	}
	}
	if( (ibx + threadIdx.x ) >= m ) 
   		tx =  -1 ;
	{
         switch(tx){
         case 0:
  	       C[0] =alpha*Cb[0] + beta * C[0];C+=ldc ; 
	 break; 
         case 1:
  	       C[0] =alpha*Cb[0] + beta * C[0];C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
	 break; 
         case 2:
  	       C[0] =alpha*Cb[0] + beta * C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
	 break; 
         case 3:
  	       C[0] =alpha*Cb[0] + beta * C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
	 break; 
         case 4:
  	       C[0] =alpha*Cb[0] + beta * C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
	 break; 
         case 5:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
	 break; 
         case 6:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
	 break; 
         case 7:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
	 break; 
         case 8:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
	 break; 
         case 9:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
	 break; 
         case 10:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
	 break; 
         case 11:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
	 break; 
         case 12:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
	 break; 
         case 13:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
	 break; 
         case 14:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
	 break; 
         case 15:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[15] + beta * C[0];C+=ldc ;
	 break; 
	default: 
		break; 
	 }
}
}
}

extern "C" __global__ void 
Ssyr2k_v16_ts_odd_generic(float *C, const float *A, const float *B, 
                          int m, int in, int k,
                          int lda, int ldb, int ldc,
                          float alpha, float beta)
{
	int tx = threadIdx.x;
	const int ty = threadIdx.y;
        int ibx = blockIdx.x ;
        int iby = blockIdx.y ;
        iby = (iby+ibx ) % gridDim.y ;
        int minor = iby &1 ;
        bool  bottom = ibx >=iby ;
        ibx = ( bottom ) ? ibx : ( iby + gridDim.y - 1 );
        iby = ( bottom ) ? iby : ( blockIdx.x + minor + gridDim.y );
        if( iby > ibx ) iby = in +1 ;
        ibx =ibx * block_M;
        iby =iby * block_N;


        const float *A1 = A ; 
        const float *B1 = B ; 

{  

	B+= iby+tx;
	B+= __mul24( ty,ldb);
	A += ibx + tx;
	C += ibx +tx +__mul24( iby+ty* unroll_f,ldc);
	float Ap[4]={A[0], A[lda], A[2*lda], A[3*lda]};

	float b=B[0];
        float b2=B[2*ldb];

	const float *Bend = B + ldb*k;
	B+=4*ldb;
	A+=4*lda;
	__shared__ float Bb[4][block_N];
	float Cb[unroll_f] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		__syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B += 4*ldb;
		__syncthreads();
	} while (B < Bend);
	Bb[ty][tx]=b;
	Bb[ty+2][tx]=b2;
	__syncthreads();
	saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);        
	saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);        
	saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);       
	saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      
	__syncthreads();

        B=A1;
        A=B1;
        int tlda = lda ; lda = ldb ; ldb = tlda ; 
	B+= iby+tx;
	B+= __mul24( ty,ldb);
	A += ibx + tx;
	Ap[0]=A[0] ;
	Ap[1]=A[lda] ;
	Ap[2]=A[2*lda] ;
	Ap[3]=A[3*lda] ;
	b=B[0];
        b2=B[2*ldb];
	const float *Bend1 = B + ldb*k;
	B+=4*ldb;
	A+=4*lda;
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		__syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B += 4*ldb;
		__syncthreads();
	} while (B < Bend1);
	Bb[ty][tx]=b;
	Bb[ty+2][tx]=b2;
	__syncthreads();
	saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);        
	saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);        
	saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);       
	saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      
	__syncthreads();

lda  = 0 ; 

if ( iby < ibx ){
tx = 15 ; 
}
else{
    if ( tx > 15 ) {
        if ( ty == 0 ){
		 lda =  1 ; 
		  tx=15;
	}
        else{ 
           lda = 1 ; 
	   tx-=16;
        }
    }
    else{
         if ( ty  == 0 ) {
		lda = 1 ; 
	}
         else {
		lda = 2 ;
		tx = 32  ;  
	}	
    }
}
if( (ibx + threadIdx.x ) >= m ) 
   tx =  -1;
{
         switch(tx){
         case 0:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
	 break; 

         case 1:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
	 break; 

         case 2:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
	 break; 

         case 3:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
	 break; 

         case 4:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
	 break; 

         case 5:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
	 break; 
         case 6:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
	 break; 
         case 7:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
	 break; 
         case 8:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
	 break; 
         case 9:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
	 break; 
         case 10:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
	 break; 
         case 11:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
	 break; 
         case 12:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
	 break; 
         case 13:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
	 break; 
         case 14:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
	 break; 
         case 15:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[15] + beta * C[0];C+=ldc ;
	 break; 
	default: 
		break; 
	 }
}

}
}

extern "C" __global__ void
Ssyr2k_v16_ts_even_special(int flag , 
                           float *C, const float *A, const float *B, 
                           int m, int in, int k, 
                           int lda, int ldb, int ldc, 
                           float alpha, float beta) 
{
	int tx = threadIdx.x;
	const int ty = threadIdx.y;
        int ibx = blockIdx.x ;
        int iby = blockIdx.y ;
        if ( flag ==1  )
        iby = (iby+ibx ) % gridDim.y ;
        const int minor = iby&1  ;
        const bool  bottom = ibx >iby ;
        ibx = ( bottom ) ? (ibx-1) : ( iby + gridDim.y  );
        iby = ( bottom ) ?  iby    : ( blockIdx.x +minor + gridDim.y );
        if( iby > ibx ) iby = in ;
        ibx =ibx * block_M;
        iby =iby * block_N;

        const float *A1 = A ; 
        const float *B1 = B ; 
{  

	B+= iby+tx;
	B+= __mul24( ty,ldb);
	A += ibx + tx;
	C += ibx +tx +__mul24( iby+ty* unroll_f,ldc);
	float Ap[4]={A[0], A[lda], A[2*lda], A[3*lda]};
	float b=B[0];
        float b2=B[2*ldb];
	const float *Bend = B + ldb*k;
	B+=4*ldb;
	A+=4*lda;
	__shared__ float Bb[4][block_N];
	float Cb[unroll_f] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		 __syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B += 4*ldb;
		 __syncthreads();
	} while (B < Bend);
	Bb[ty][tx]=b;
	Bb[ty+2][tx]=b2;
	__syncthreads();
	saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);        
	saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);        
	saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);       
	saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      

// -- 2nd Half            
        B=A1;
        A=B1;
        int tlda = lda ; lda = ldb ; ldb = tlda ; 

	B+= iby+tx;
	B+= __mul24( ty,ldb);
	A += ibx + tx;
	Ap[0]=A[0] ;
	Ap[1]=A[lda] ;
	Ap[2]=A[2*lda] ;
	Ap[3]=A[3*lda] ;
	b=B[0];
        b2=B[2*ldb];
	const float *Bend1 = B + ldb*k;
	B+=4*ldb;
	A+=4*lda;
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		 __syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B += 4*ldb;
		 __syncthreads();
	} while (B < Bend1);
	Bb[ty][tx]=b;
	Bb[ty+2][tx]=b2;
	__syncthreads();
	saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);        
	saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);        
	saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);       
	saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      

lda  = 0 ; 

if ( iby < ibx ){
/*	#pragma unroll 16 
         for (int i = 0; i < unroll_f; i++, C += ldc) 
  	       C[0] =alpha*Cb[i] + beta * C[0];
*/
tx = 15 ; 
}
else{
    if ( tx > 15 ) {
        if ( ty == 0 ){
		 lda =  1 ; 
		  tx=15;
	}
        else{ 
           lda = 1 ; 
	   tx-=16;
        }
    }
    else{
         if ( ty  == 0 ) {
		lda = 1 ; 
	}
         else {
		lda = 2 ;
		tx = 32  ;  
	}	
    }
}
{
         switch(tx){
         case 0:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
	 break; 
         case 1:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
	 break; 
         case 2:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
	 break; 
         case 3:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
	 break; 
         case 4:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
	 break; 
         case 5:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
	 break; 
         case 6:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
	 break; 
         case 7:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
	 break; 
         case 8:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
	 break; 
         case 9:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
	 break; 
         case 10:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
	 break; 
         case 11:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
	 break; 
         case 12:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
	 break; 
         case 13:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
	 break; 
         case 14:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
	 break; 
         case 15:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[15] + beta * C[0];C+=ldc ;
	 break; 
	default: 
		break; 
	 }
}

}
}

extern "C" __global__ void 
Ssyr2k_v16_ts_odd_special(int flag,
                          float *C, const float *A, const float *B, 
                          int m, int in, int k,
                          int lda, int ldb, int ldc,
                          float alpha, float beta) 
{
	int tx = threadIdx.x;
	const int ty = threadIdx.y;
        int ibx = blockIdx.x ;
        int iby = blockIdx.y ;
        if ( flag ==1  )
        iby = (iby+ibx ) % gridDim.y ;
        int minor = iby &1 ;
        bool  bottom = ibx >=iby ;
        ibx = ( bottom ) ? ibx : ( iby + gridDim.y - 1 );
        iby = ( bottom ) ? iby : ( blockIdx.x + minor + gridDim.y );
        if( iby > ibx ) iby = in +1 ;
        ibx =ibx * block_M;
        iby =iby * block_N;


        const float *A1 = A ; 
        const float *B1 = B ; 

if( iby > ibx) {
 return ; 
}
else{  

	B+= iby+tx;
	B+= __mul24( ty,ldb);
	A += ibx + tx;
	C += ibx +tx +__mul24( iby+ty* unroll_f,ldc);
	float Ap[4]={A[0], A[lda], A[2*lda], A[3*lda]};

	float b=B[0];
        float b2=B[2*ldb];

	const float *Bend = B + ldb*k;
	B+=4*ldb;
	A+=4*lda;
	__shared__ float Bb[4][block_N];
	float Cb[unroll_f] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		__syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B += 4*ldb;
		__syncthreads();
	} while (B < Bend);
	Bb[ty][tx]=b;
	Bb[ty+2][tx]=b2;
	__syncthreads();
	saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);        
	saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);        
	saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);       
	saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      
        B=A1;
        A=B1;
        int tlda = lda ; lda = ldb ; ldb = tlda ; 
	B+= iby+tx;
	B+= __mul24( ty,ldb);
	A += ibx + tx;
	Ap[0]=A[0] ;
	Ap[1]=A[lda] ;
	Ap[2]=A[2*lda] ;
	Ap[3]=A[3*lda] ;
	b=B[0];
        b2=B[2*ldb];
	const float *Bend1 = B + ldb*k;
	B+=4*ldb;
	A+=4*lda;
	do {
		float Ab[4] = {Ap[0], Ap[1], Ap[2], Ap[3]};
		Bb[ty][tx]=b;
		Bb[ty+2][tx]=b2;
		__syncthreads();
		Ap[0] = A[0];
		Ap[1] = A[lda];
		Ap[2] = A[2*lda];
		Ap[3] = A[3*lda];
		b=B[0];
   	        b2=B[2*ldb];
		saxpy(Ab[0], &Bb[0][ty*unroll_f], Cb);
		saxpy(Ab[1], &Bb[1][ty*unroll_f], Cb);
		saxpy(Ab[2], &Bb[2][ty*unroll_f], Cb);
		saxpy(Ab[3], &Bb[3][ty*unroll_f], Cb);
		A+=4*lda;
		B += 4*ldb;
		__syncthreads();
	} while (B < Bend1);
	Bb[ty][tx]=b;
	Bb[ty+2][tx]=b2;
	__syncthreads();
	saxpy(Ap[0], &Bb[0][ty*unroll_f], Cb);        
	saxpy(Ap[1], &Bb[1][ty*unroll_f], Cb);        
	saxpy(Ap[2], &Bb[2][ty*unroll_f], Cb);       
	saxpy(Ap[3], &Bb[3][ty*unroll_f], Cb);      

lda  = 0 ; 

if ( iby < ibx ){
/*	#pragma unroll 16 
         for (int i = 0; i < unroll_f; i++, C += ldc) 
  	       C[0] =alpha*Cb[i] + beta * C[0];
*/
tx = 15 ; 
}
else{
    if ( tx > 15 ) {
        if ( ty == 0 ){
		 lda =  1 ; 
		  tx=15;
	}
        else{ 
           lda = 1 ; 
	   tx-=16;
        }
    }
    else{
         if ( ty  == 0 ) {
		lda = 1 ; 
	}
         else {
		lda = 2 ;
		tx = 32  ;  
	}	
    }
}

{
         switch(tx){
         case 0:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
	 break; 

         case 1:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
	 break; 

         case 2:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
	 break; 

         case 3:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
	 break; 

         case 4:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
	 break; 

         case 5:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
	 break; 
         case 6:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
	 break; 
         case 7:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
	 break; 
         case 8:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
	 break; 
         case 9:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
	 break; 
         case 10:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
	 break; 
         case 11:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
	 break; 
         case 12:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
	 break; 
         case 13:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
	 break; 
         case 14:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
	 break; 
         case 15:
  	       C[0] =alpha*Cb[0] + beta*C[0]; C+=ldc ; 
  	       C[0] =alpha*Cb[1] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[2] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[3] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[4] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[5] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[6] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[7] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[8] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[9] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[10] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[11] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[12] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[13] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[14] + beta * C[0];C+=ldc ;
  	       C[0] =alpha*Cb[15] + beta * C[0];C+=ldc ;
	 break; 
	default: 
		break; 
	 }
}

}
}

#include<stdio.h>

extern "C" void
magmablas_ssyr2k(char UPLO, char TRANS, int m , int k, float alpha,  
                 const float *A, int lda , const float *B, int ldb, 
                 float beta, float *C, int ldc)
{
/* -- MAGMA (version 1.0) --
      Univ. of Tennessee, Knoxville
      Univ. of California, Berkeley
      Univ. of Colorado, Denver
      November 2010

   Purpose
   =======
   SSYR2K  performs one of the symmetric rank 2k operations
      C := alpha*A*B' + alpha*B*A' + beta*C,
   or
      C := alpha*A'*B + alpha*B'*A + beta*C,

   where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
   and  A and B  are  n by k  matrices  in the  first  case  and  k by n
   matrices in the second case.

   This implementation is for UPLO == 'L' and TRANS == 'N'.
  
   Assumptions
   ===========
   Both lda and ldb must be multiple of 32. 
   Parameter k must be divisible by 8 - note that this algorithm was developed
   for the tridiagonal factorization and k in that case would be the blocking size.
   We always request the blocking size to be divisible by at least 16.

   This kernel goes to about 300 GFlop/s on the GTX280.
   ====================================================================== */	

    int in = m / block_M;
    int flag = 1 ;
    if ( lda >=1024 && lda %256 == 0 )
      flag = 1 ; // It was kept to reorder the GPUs internal scheduling of thread blocks. 
      if( m  % block_M == 0 ) {
        if ( in&1 )
          {
            dim3 grid( in, (in/2+1));
            dim3 threads( thread_x, thread_y );
            Ssyr2k_v16_ts_odd_special<<< grid, threads >>>(flag,  
                                                           C, A, B, 
                                                           m, in/2, k,
                                                           lda, ldb, ldc, 
                                                           alpha, beta);
         }
       else
        {
           dim3 grid( in+1, (in/2));
           dim3 threads( thread_x, thread_y );
           Ssyr2k_v16_ts_even_special<<< grid, threads >>>(flag,  
                                                           C, A, B, 
                                                           m, in/2, k, 
                                                           lda, ldb, ldc,
                                                           alpha, beta);
       }
     }
   else{
    in+=1;
    if( in&1 )
      {  
         dim3 grid( in, (in/2+1));
         dim3 threads( thread_x, thread_y );
         Ssyr2k_v16_ts_odd_generic<<< grid, threads >>>(C, A, B,
                                                        m, in/2, k, 
                                                        lda, ldb, ldc,
                                                        alpha, beta);
      }
    else
      {
         dim3 grid( in+1, (in/2));
         dim3 threads( thread_x, thread_y );
         Ssyr2k_v16_ts_even_generic<<< grid, threads >>>(C, A, B,
                                                         m, in/2, k, 
                                                         lda, ldb, ldc,
                                                         alpha, beta);
      }
   }
}

