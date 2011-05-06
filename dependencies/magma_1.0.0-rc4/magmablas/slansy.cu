/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated s

*/
#include "common_magma.h"

#define BLOCK_SIZE 32

//#define num_threads 64
#define dgemv_bs 32

#define PRECISION_s
#if (!defined(PRECISION_z)) || (GPUSHMEM >= 200)

__global__ void
l_slansy_special (int n, float* A, int lda,  float *y){
  int tx = threadIdx.x ; 
  int ty = threadIdx.y ; 
  int ind = blockIdx.x*  dgemv_bs + tx ;
  float res = 0.;

  __shared__ float la[dgemv_bs][dgemv_bs+1];
  	
  A += ind;
  A+= ty * lda  ;  
  int break_d  =   blockIdx.x* dgemv_bs ;

  for(int  i=0; i<break_d; i += dgemv_bs ){
    #pragma unroll 8 
    for(int j=0; j < dgemv_bs ; j+=4){
        la[tx][ty+j] = A[j*lda] ;
    }
    __syncthreads();

    #pragma unroll 8 
    for(int j=0; j < 8 ; j++){
       res+=fabsf( la[tx][j+ty*8]) ;
    }
    A+=lda* dgemv_bs ;
    __syncthreads(); 
  }

 
  #pragma unroll 8
  for(int j =0; j<dgemv_bs; j+=4)
         la[ty+j][tx] = A[ j * lda];


  A+= dgemv_bs ;
  __syncthreads();
  #pragma unroll 8
  for(int  i=ty*8; i<(1+ty)* dgemv_bs/4 ; i++){
         if ( i < tx )   {
	        la[tx][i] = la[i][tx] ; 
         }
	 else 
	        la[tx][i] = la[tx][i]  ;
  
  }
  __syncthreads();
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4 ; j++){
     res+=fabsf(la[tx][j+ty*8]);
    }
  break_d  += dgemv_bs ; 
  __syncthreads();

  for(int i=break_d; i<n; i += dgemv_bs ){
   #pragma unroll 8
    for(int j=0; j<dgemv_bs; j+=4)
       la[ty+j][tx] = A[ j * lda];
    A+= dgemv_bs ;
      __syncthreads();
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4;j++){
       res+= fabsf(la[tx][j+ty*8]);
    }
      __syncthreads();
  }


  la[tx][ty] = MAGMA_S_MAKE( res, 0. );
   __syncthreads();
   if( ty == 0 ) {
     res = res 
       + MAGMA_S_REAL( la[tx][1] ) 
       + MAGMA_S_REAL( la[tx][2] )
       + MAGMA_S_REAL( la[tx][3] );
     y[ind] = res;
   }

}

__global__ void
l_slansy_generic(int n, float* A, int lda,  float *y , int m_full_block , 
                 int m_mod_32)
{ 
  int tx = threadIdx.x ; 
  int ty = threadIdx.y ; 

  int ind = blockIdx.x*  dgemv_bs + tx ;
  
  float res = 0.;

  __shared__ float la   [dgemv_bs][dgemv_bs+1];

  if( blockIdx.x == m_full_block ) {
  /************************************************************************
   -- Last block --
   -- We will do something unusual here 
   -- For sufficiently large matrix the overhead will be very low
  *************************************************************************/
       if  ( tx < m_mod_32 ){
		A+= ( blockIdx.x * dgemv_bs + tx ) ;
       } 	 	
       else{
		A+= ( blockIdx.x * dgemv_bs + m_mod_32 -1) ; 
       }
       A+= ty * lda  ;  
       int break_d  =   blockIdx.x* dgemv_bs ;

	  /*----------------------------
		Go Right
	  -------------------------------*/

	  for(int  i=0; i<break_d; i += dgemv_bs ){
	    #pragma unroll 8 
	    for(int j=0; j < dgemv_bs ; j+=4){
	        la[tx][ty+j] = A[j*lda] ;
	    }
	    __syncthreads();

	    #pragma unroll 8 
	    for(int j=0; j < 8 ; j++){
	       res+=fabsf( la[tx][j+ty*8]);
	    }
	    A+=lda* dgemv_bs ;
	    __syncthreads(); 
	  }
	  /*
           we don't need to make zero, as those computation will be discarded. 
          */
          if( ty==0  ) {
		/*--------------------------------------------
			he will compute the triangular parts
			others will be waiting with values. 
                -----------------------------------------------*/
		int j ;
                int count = 1 ; 
		if( tx < m_mod_32 ) 
			count = tx ; 
		else
			count = m_mod_32 ;
		for(j =0;j<=count;j++){
			res+= fabsf( A[j*lda]) ;
                }
		A+=(tx)*lda;
		count = 1 ; 
		for(;j<m_mod_32;j++){
			res+=fabsf( A[count]) ;
			count++;
		}
          }
          else{
          }
	  __syncthreads(); 
          la[tx][ty]= MAGMA_S_MAKE( res, 0. ) ;
          __syncthreads();
         /*--------------------------------------------------------
	 The leader accumulates all the results from his peer. 
         ----------------------------------------------------------*/
         if( ty == 0 ) {
           res = res 
             + MAGMA_S_REAL( la[tx][1] ) 
             + MAGMA_S_REAL( la[tx][2] )
             + MAGMA_S_REAL( la[tx][3] );
           if( tx < m_mod_32)
             y[ind] = res;
         }
	 
  }

  else{ 
  /***************************************
    -----------------------------------
  -- All the blocks but the last one --
  ****************************************
  -------------------------------------*/
  A += ind;
  A+= ty * lda  ;  
  int break_d  =   blockIdx.x* dgemv_bs ;

  /*----------------------------
	Go Right
  -------------------------------*/
  for(int  i=0; i<break_d; i += dgemv_bs ){
    #pragma unroll 8 
    for(int j=0; j < dgemv_bs ; j+=4){
        la[tx][ty+j] = A[j*lda] ;
    }
    __syncthreads();

    #pragma unroll 8 
    for(int j=0; j < 8 ; j++){
       res+=fabsf(la[tx][j+ty*8]);
    }
    A+=lda* dgemv_bs ;
    __syncthreads(); 
  }

 
  /*------------------------------------
	Diagonal 
	Copy + Transpose lower triangle
  --------------------------------------*/
  #pragma unroll 8
  for(int j =0; j<dgemv_bs; j+=4)
         la[ty+j][tx] = A[ j * lda];


  A+= dgemv_bs ;
  __syncthreads();
  /*--------------------------------------------
	Mirror Upper Triangle to Lower triangle
  ---------------------------------------------*/
  #pragma unroll 8
  for(int  i=ty*8; i<(1+ty)* dgemv_bs/4 ; i++){
         if ( i < tx )   {
	        la[tx][i] = la[i][tx] ; 
         }
	 else 
	        la[tx][i] = la[tx][i]  ;
  
  }
  __syncthreads();
  /*--------------------------------
	Do diagonal Computation
  -----------------------------------*/
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4 ; j++){
     res+= fabsf(la[tx][j+ty*8]);
    }
  break_d  += dgemv_bs ; 
  __syncthreads();


  n -= m_mod_32 ;  // @ 
  /*-----------------------------
	Go Down 
  -------------------------------*/
  for(int i=break_d; i<n; i += dgemv_bs ){
   #pragma unroll 8
    for(int j=0; j<dgemv_bs; j+=4)
       la[ty+j][tx] = A[ j * lda];
    A+= dgemv_bs ;
      __syncthreads();
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4;j++){
       res+=fabsf(la[tx][j+ty*8]);
    }
      __syncthreads();
  }

  
  /*---------------------------------------------
	doing m_mod_32 stuffs here.
	Symmetric is giving us benefit .. true
  -----------------------------------------------*/
    A-=tx;
    if( tx < m_mod_32){
	A+=tx;
    }
    else{
	A+=(m_mod_32-1); /* Same as above*/
    }

   #pragma unroll 8
    for(int j=0; j<dgemv_bs; j+=4){
       if( tx < m_mod_32 ) 
         la[ty+j][tx] = MAGMA_S_MUL( MAGMA_S_ONE,  A[ j * lda] );
       else
         la[ty+j][tx] = MAGMA_S_MUL( MAGMA_S_ZERO, A[ j * lda] );
       
    }
    __syncthreads();

    /*----------------------------------------
	What about doing some Zeroing here?
	instead of zeroing before?
    -----------------------------------------*/	
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4;j++){
       res+=fabsf(la[tx][j+ty*8]);
    }
    __syncthreads();
   

    la[tx][ty]= MAGMA_S_MAKE( res, 0. );
   __syncthreads();
   /*--------------------------------------------------------
	The leader accumulates all the results from his peer. 
   ----------------------------------------------------------*/
   if( ty == 0 ) {
     res = res 
       + MAGMA_S_REAL( la[tx][1] ) 
       + MAGMA_S_REAL( la[tx][2] )
       + MAGMA_S_REAL( la[tx][3] );
     y[ind] = res;
   }

  }

}

__global__ void
u_slansy_generic (int n, float* A, int lda, float *y , int m_full_block , int m_mod_32){

  
  int tx = threadIdx.x ; 
  int ty = threadIdx.y ; 

  int ind = blockIdx.x*  dgemv_bs + tx ;
  
  float res = 0.;


  __shared__ float la   [dgemv_bs][dgemv_bs+1];
  int blockIdxx =  blockIdx.x ;

  if( blockIdx.x == m_full_block ) {

  /************************************************************************
   -- Last block --
   -- We will do something unusual here 
   -- For sufficiently large matrix the overhead will be very low
  *************************************************************************/

  ind =  tx ;
  A+= lda*(n-1) ; 


       if  ( tx < m_mod_32 ){
		A+= (  tx ) ;
       } 	 	
       else{
		A+= (  m_mod_32 -1) ; 
       }
       A-= ty * lda  ;  
       int break_d  =   (blockIdx.x)* dgemv_bs ;

	  /*----------------------------
		Go Right
	  -------------------------------*/

	  for(int  i=0; i<break_d; i += dgemv_bs ){
	    #pragma unroll 8 
	    for(int j=0; j < dgemv_bs ; j+=4){
	        la[tx][ty+j] = A[-j*lda] ;
	    }
	    __syncthreads();

	    #pragma unroll 8 
	    for(int j=0; j < 8 ; j++){
	       res+=fabsf(la[tx][j+ty*8]);
	    }
	    A-=lda* dgemv_bs ;
	    __syncthreads(); 
	  }
	  /*
           we don't need to make zero, as those computation will be discarded. 
          */
          if( ty==0  ) {
		/*--------------------------------------------
			he will compute the triangular parts
			others will be waiting with values. 
                -----------------------------------------------*/
		int j ;
                int count = 1 ; 
		if( tx < m_mod_32 ) 
			count =m_mod_32- tx ; 
		else
			count = m_mod_32 ;
		for(j =0;j<count;j++){
			res+= fabsf( A[-j*lda] );
                }
		A-=(count-1)*lda;
		count = 1 ; 
		for(;j<m_mod_32;j++){
			res+= fabsf( A[-count] );
			count++;
		}
          }
          else{
          }
	  __syncthreads(); 
          la[tx][ty]= MAGMA_S_MAKE( res, 0. );
          __syncthreads();
         /*--------------------------------------------------------
	 The leader accumulates all the results from his peer. 
         ----------------------------------------------------------*/
         if( ty == 0 ) {
           res = res 
             + MAGMA_S_REAL( la[tx][1] ) 
             + MAGMA_S_REAL( la[tx][2] )
             + MAGMA_S_REAL( la[tx][3] );
           if( tx < m_mod_32)
             y[ind] = res;
         }
	 
  }

  else{ 
  /***************************************
    -----------------------------------
  -- All the blocks but the last one --
  -- By the way this code can be optimized more. 
  ****************************************
  -------------------------------------*/
  ind = blockIdx.x *  dgemv_bs + tx + m_mod_32 ;
  float *A1 = A ; 
  A+= lda*(n-1)  ; 

  A += ind;
  A-= ty * lda  ;  

  int break_d  = (n / dgemv_bs -   blockIdxx-1 )* dgemv_bs ;
  /*----------------------------
	Go Left
  -------------------------------*/
  for(int  i=0; i<break_d; i += dgemv_bs ){
    #pragma unroll 8 
    for(int j=0; j < dgemv_bs ; j+=4){
        la[tx][ty+j] = A[-j*lda] ;
    }
    __syncthreads();

    #pragma unroll 8 
    for(int j=0; j < 8 ; j++){
       res+=fabsf( la[tx][j+ty*8]);
    }
    A-=lda* dgemv_bs ;
    __syncthreads(); 
  }

 
  /*------------------------------------
	Diagonal 
	Copy + Transpose lower triangle
  --------------------------------------*/
  #pragma unroll 8
  for(int j =0; j<dgemv_bs; j+=4){
         la[tx][31-ty-j] = A[ -j * lda];
  }

  A-= dgemv_bs ;
  __syncthreads();
  /*--------------------------------------------
	Mirror Upper Triangle to Lower triangle
  ---------------------------------------------*/
  #pragma unroll 8
  for(int  i=ty*8; i<(1+ty)* dgemv_bs/4 ; i++){
         if ( i <tx ){
	        la[tx][i] = la[i][tx]; 
         }
	 else{ 
	        la[tx][i] = la[tx][i]  ;
	 }
  }
  __syncthreads();
  /*--------------------------------
	Do diagonal Computation
  -----------------------------------*/
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4 ; j++){
     res+=fabsf(  la[tx][j+ty*8] ) ;
    }
  break_d  += dgemv_bs ; 
  __syncthreads();


  n -= m_mod_32 ;  // @ 
  /*-----------------------------
	Go Up 
  -------------------------------*/
  int i ;
  for( i=break_d; i<n; i+= dgemv_bs ){
   #pragma unroll 8
    for(int j=0; j<dgemv_bs; j+=4){
       la[ty+j][tx] = A[- j * lda];
    }
    A-= dgemv_bs ;
      __syncthreads();
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4;j++){
       res+=fabsf ( la[31-tx][j+ty*8] );
    }
      __syncthreads();
  }
  /*---------------------------------------------
	doing m_mod_32 stuffs here.
	Symmetric is giving us benefit .. true
	Do the other way please......
  -----------------------------------------------*/
   A1 = A1 + m_mod_32 * lda + tx *lda ;  
   if( ty == 0  ) {
	for( int j = 0 ;  j < m_mod_32 ; j++){
		res+=  fabsf (  A1[ j + lda * (blockIdx.x) * 32 ] ) ;
	}
   }
    __syncthreads();

    la[tx][ty]= MAGMA_S_MAKE( res, 0);
   __syncthreads();
   /*--------------------------------------------------------
	The leader accumulates all the results from his peer. 
   ----------------------------------------------------------*/
   if( ty == 0 ) {
     res = res 
       + MAGMA_S_REAL( la[tx][1] ) 
       + MAGMA_S_REAL( la[tx][2] )
       + MAGMA_S_REAL( la[tx][3] );
     y[ind] =  res;
   }

  }

}

__global__ void
u_slansy_special (int n, float* A, int lda, float *y ){
  int tx = threadIdx.x ; 
  int ty = threadIdx.y ; 
  int ind = blockIdx.x*  dgemv_bs + tx ;
  float res = 0.;

  /*
	Reverse Computation ... 
		- Left 
		- Triangle 
		- Up 
  */

  A+= lda*(n-1) ; 
  __shared__ float la   [dgemv_bs][dgemv_bs+1];


  A += ind;
  A-= ty * lda  ;  
  int break_d  = (n / dgemv_bs -   blockIdx.x-1 )* dgemv_bs ;

  for(int  i=0; i<break_d; i += dgemv_bs ){
    #pragma unroll 8 
    for(int j=0; j < dgemv_bs ; j+=4){
        la[tx][ty+j] = A[-j*lda] ;
    }
    __syncthreads();

    #pragma unroll 8 
    for(int j=0; j < 8 ; j++){
       res+=fabsf(la[tx][j+ty*8]);
    }
    A-=lda* dgemv_bs ;
    __syncthreads(); 
  }




  #pragma unroll 8
  for(int j =0; j<dgemv_bs; j+=4)
         la[tx][31-ty-j] = A[ -j * lda];
  /*
	Look at the indexing changes
  */

  A-= dgemv_bs ;
  __syncthreads();
  #pragma unroll 8
  for(int  i=ty*8; i<(1+ty)* dgemv_bs/4 ; i++){
         if ( i <tx ){
	        la[tx][i] = la[i][tx]; 
         }
	 else{ 
	        la[tx][i] = la[tx][i]  ;
	 }
  
  }
  __syncthreads();
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4 ; j++){
     res+= fabsf(la[tx][j+ty*8]);
    }

  break_d  += dgemv_bs ; 
  __syncthreads();



  for(int i=break_d; i<n; i+= dgemv_bs ){
   #pragma unroll 8
    for(int j=0; j<dgemv_bs; j+=4)
       la[ty+j][tx] = A[ -j * lda];

    A-= dgemv_bs ;
      __syncthreads();
    #pragma unroll 8
    for(int j=0; j < dgemv_bs/4;j++){
       res+=fabsf( la[31-tx][j+ty*8]);
    }
      __syncthreads();
  }


  la[tx][ty]= MAGMA_S_MAKE( res, 0. );

   __syncthreads();
   if( ty == 0 ) {
     res = res 
       + MAGMA_S_REAL( la[tx][1] ) 
       + MAGMA_S_REAL( la[tx][2] )
       + MAGMA_S_REAL( la[tx][3] );
     y[ind] =   res;
   }

}





extern "C" void mslansy (char uplo , int m ,  float *A , int lda ,  float *Y  )
{
/*
Note:
	The UPLO = 'U' Version can be optimized more.
        side is not needed........................... 
*/
    int blocks;
    if (m % dgemv_bs==0)
        blocks = m/ dgemv_bs;
    else
        blocks = m/ dgemv_bs + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(32, 4, 1);

    if( m % dgemv_bs == 0 ) {
	    if( uplo == 'L' || uplo == 'l'){	
		    l_slansy_special <<<grid, threads>>> (m, A, lda, Y);
	    }
            else{
		    u_slansy_special <<<grid, threads>>> (m, A, lda,  Y);
	    } 
		
    } 
    else{	
	    int  m_full_block = (m - m % 32 ) /32 ; 
	    int  m_mod_32 = m%32 ;  
	    if( uplo == 'L' || uplo == 'l'){
		    l_slansy_generic <<<grid, threads>>> (m, A, lda, Y , m_full_block , m_mod_32);
	    }	
	    else{
		    u_slansy_generic <<<grid, threads>>> (m, A, lda, Y , m_full_block , m_mod_32);
	    }	
    }
}


extern "C" float 
magmablas_slansy(char norm, char uplo, int n, 
                 float *A, int lda, float *WORK )
{
	if( norm != 'I' && norm !='i')	{
                printf("Only normI is available\n");
		exit(1);
	}
        else{
		mslansy ( uplo , n , A , lda , WORK);
		int val = cublasIsamax(n,WORK,1);
	        float retVal[1];
		cublasGetMatrix( 1, 1, sizeof( float ), WORK+val-1, 1, retVal, 1 ) ;
		return retVal[0];
        }

}

#endif /* (!defined(PRECISION_z)) || (GPUSHMEM >= 200) */
