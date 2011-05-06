/*
  -- MAGMA (version 1.0) --
  Univ. of Tennessee, Knoxville
  Univ. of California, Berkeley
  Univ. of Colorado, Denver
  November 2010

  @precisions mixed zc -> ds

*/
#include "common_magma.h"
#define PRECISION_z

#define BLOCK_SIZE 32
//#define num_threads 64
#define dgemv_bs 32

#if (!defined(PRECISION_z)) || (GPUSHMEM >= 200)

/*------------------------------------------ UPLO = 'L' ----------------------------------*/
__device__ int flag = 0; 

__global__ void
l_zlat2c_special (int n, cuDoubleComplex* A, int lda,  cuFloatComplex *SA , int *INFO, double RMAX, int ldsa )
{
    double mRMAX = - RMAX;
    int tx  = threadIdx.x ;
    int ty  = threadIdx.y ;
    int ind = blockIdx.x*  dgemv_bs + tx ;

    __shared__ cuDoubleComplex la[dgemv_bs][dgemv_bs+1];

    A  += ind;
    SA += ind;
    A  += ty * lda;
    SA += ty * ldsa;

    int break_d  =   blockIdx.x* dgemv_bs ;
    cuDoubleComplex temp ;
    for(int  i=0; i<break_d; i += dgemv_bs ){
        #pragma unroll 8
        for(int j=0; j < dgemv_bs; j+=4){
            temp = A[j*lda] ;
            if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                )
                {
                    flag = 1;
                }
            SA[j*ldsa] = cuComplexDoubleToFloat(temp);
        }
        A  += lda *dgemv_bs ;
        SA += ldsa*dgemv_bs ;
    }

#pragma unroll 8
    for(int j =0; j<dgemv_bs; j+=4)
        la[ty+j][tx] = A[j*lda];

    __syncthreads();

    A += dgemv_bs ;
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
    for(int j =0; j<dgemv_bs; j+=4){
        temp =  la[ty+j][tx] ;
        if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
            || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
            )
            {
                flag = 1;
            }
        SA[j*ldsa] = cuComplexDoubleToFloat(temp);
    }

    __syncthreads();
    //la[tx][ty] = flag ;
    //__syncthreads();
    //if( ty == 0 ) {
    //    //  INFO[0] = flag+ la[tx] [1] +  la[tx] [2] + la[tx] [3] ;
    //}
}

__global__ void
l_zlat2c_generic(int n, cuDoubleComplex* A, int lda, cuFloatComplex *SA , int m_full_block,
                 int m_mod_32 , int *INFO , double RMAX , int ldsa)
{
    double mRMAX = - RMAX;
    int tx = threadIdx.x ;
    int ty = threadIdx.y ;

    int ind = blockIdx.x*  dgemv_bs + tx ;


    __shared__ cuDoubleComplex la   [dgemv_bs][dgemv_bs+1];

    cuDoubleComplex temp ;

    if( blockIdx.x == m_full_block ) {
        /************************************************************************
   -- Last block --
   -- We will do something unusual here
   -- For sufficiently large matrix the overhead will be very low
        *************************************************************************/
        if  ( tx < m_mod_32 ){
            A+= ( blockIdx.x * dgemv_bs + tx ) ;
            SA+= ( blockIdx.x * dgemv_bs + tx ) ;
        }
        else{
            A+= ( blockIdx.x * dgemv_bs + m_mod_32 -1) ;
            SA+= ( blockIdx.x * dgemv_bs + m_mod_32 -1) ;
        }
        A+= ty * lda  ;
        SA+= ty * ldsa  ;
        int break_d  =   blockIdx.x* dgemv_bs ;

        /*----------------------------
          Go Right
          -------------------------------*/

        for(int  i=0; i<break_d; i += dgemv_bs ){
#pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4){
                temp = A[j*lda] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[j*ldsa] = cuComplexDoubleToFloat(temp);
            }
            A  += lda *dgemv_bs ;
            SA += ldsa*dgemv_bs ;
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
                temp =  A[j*lda] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[j*ldsa] = cuComplexDoubleToFloat(temp);
            }
            A  += (tx)*lda;
            SA += (tx)*ldsa;
            count = 1 ;
            for(;j<m_mod_32;j++){
                temp= A[count] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[count] = cuComplexDoubleToFloat(temp);
                count++;
            }
        }
        else{
        }

        //la[tx][ty] = flag ;

        __syncthreads();
        /*--------------------------------------------------------
          The leader accumulates all the results from his peer.
          ----------------------------------------------------------*/
        //if( ty == 0 ) {
            //  INFO[ind] = ld[tx][0] +  ld[tx][1] + ld[tx][2] + ld[tx][3] ;
        //}

    }

    else{
        /***************************************
    -----------------------------------
  -- All the blocks but the last one --
  ****************************************
  -------------------------------------*/
        A += ind;
        SA += ind;
        A+= ty * lda  ;
        SA+= ty * ldsa  ;
        int break_d  =   blockIdx.x* dgemv_bs ;

        /*----------------------------
          Go Right
          -------------------------------*/
        for(int  i=0; i<break_d; i += dgemv_bs ){
#pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4){
                temp = A[j*lda] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[j*ldsa] = cuComplexDoubleToFloat(temp);
            }
            A  += lda *dgemv_bs ;
            SA += ldsa*dgemv_bs ;
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
        for(int j =0; j<dgemv_bs; j+=4){
            temp =  la[ty+j][tx] ;
            if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                )
                {
                    flag = 1;
                }
            SA[j*ldsa] = cuComplexDoubleToFloat(temp);
        }
        __syncthreads();
        //la[tx] [ty ] = flag ;
        //__syncthreads();
        /*--------------------------------------------------------
          The leader accumulates all the results from his peer.
          ----------------------------------------------------------*/
        //if( ty == 0 )
        //    {
                //  INFO [ind] =  flag + la[tx][1]+ la[tx][2]+ la[tx][3] ;
        //    }
    }
}


/*- ---------------------------------------------- UPLO = 'U' ----------------------------------*/
/* Generic Case*/

__global__ void
u_zlat2c_generic(int n, cuDoubleComplex* A, int lda, cuFloatComplex *SA, int m_full_block,
                 int m_mod_32 , int *INFO , double RMAX , int ldsa)
{
    double mRMAX = - RMAX;
    int tx = threadIdx.x ;
    int ty = threadIdx.y ;

    int ind = blockIdx.x*  dgemv_bs + tx ;


    __shared__ cuDoubleComplex la   [dgemv_bs][dgemv_bs+1];
    int blockIdxx =  blockIdx.x ;
    cuDoubleComplex temp ;
    if( blockIdx.x == m_full_block ) {

        /************************************************************************
   -- Last block --
   -- We will do something unusual here
   -- For sufficiently large matrix the overhead will be very low
        *************************************************************************/

        ind =  tx ;
        A+= lda*(n-1) ;
        SA+= ldsa*(n-1) ;


        if  ( tx < m_mod_32 ){
            A+= (  tx ) ;
            SA+= (  tx ) ;
        }
        else{
            A+= (  m_mod_32 -1) ;
            SA+= (  m_mod_32 -1) ;
        }
        A-= ty * lda  ;
        SA-= ty * ldsa  ;
        int break_d  =   (blockIdx.x)* dgemv_bs ;

        /*----------------------------
          Go Right
          -------------------------------*/

        for(int  i=0; i<break_d; i += dgemv_bs ){
#pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4){
                temp  = A[-j*lda] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[-j*ldsa] = cuComplexDoubleToFloat(temp);
            }
            A  -= lda *dgemv_bs ;
            SA -= ldsa*dgemv_bs ;
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
                temp =  A[-j*lda] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[-j*ldsa] = cuComplexDoubleToFloat(temp);
            }
            A-=(count-1)*lda;
            SA-=(count-1)*ldsa;
            count = 1 ;
            for(;j<m_mod_32;j++){
                temp =  A[-count] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[-count] = cuComplexDoubleToFloat(temp);
                count++;
            }
        }
        else{
        }

        /*--------------------------------------------------------
          The leader accumulates all the results from his peer.
          ----------------------------------------------------------*/
        //la[tx][ty] = flag ;
        //__syncthreads();
        //if( ty == 0 ) {
        //    // INFO [ind] = flag  + la[tx][1] + la[tx][2] + la[tx][3]  ;
        //}

    }

    else{
        /***************************************
    -----------------------------------
  -- All the blocks but the last one --
  -- By the way this code can be optimized more.
  ****************************************
  -------------------------------------*/
        ind = blockIdx.x *  dgemv_bs + tx + m_mod_32 ;
        A+= lda*(n-1)  ;
        SA+= ldsa*(n-1)  ;


        A += ind;
        SA += ind;
        A-= ty * lda  ;
        SA-= ty * ldsa  ;

        int break_d  = (n / dgemv_bs -   blockIdxx-1 )* dgemv_bs ;
        /*----------------------------
          Go Left
          -------------------------------*/
        for(int  i=0; i<break_d; i += dgemv_bs ){
#pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4){
                temp  = A[-j*lda] ;
                if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                    || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                    )
                    {
                        flag = 1;
                    }
                SA[-j*ldsa] = cuComplexDoubleToFloat(temp);
            }
            A-=lda* dgemv_bs ;
            SA-=ldsa* dgemv_bs ;
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
        for(int j =0; j<dgemv_bs; j+=4){
            temp =   la[tx][31-ty-j];
            // temp =  la[ty+j][tx] ;
            if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                )
                {
                    flag = 1;
                }
            SA[- j*ldsa] = cuComplexDoubleToFloat(temp);
        }

        /*--------------------------------------------------------
          The leader accumulates all the results from his peer.
          ----------------------------------------------------------*/
        //la[tx] [ty] = flag ;
        //__syncthreads();
        //
        //if( ty == 0 ) {
        //    //  INFO[ind] = flag +  la[tx] [1] +  la[tx] [2] + la[tx] [3]  ;
        //}

    }
}

/*- ---------------------------------------------- UPLO = 'U' ----------------------------------*/
/*Good Dimension*/
__global__ void
u_zlat2c_special (int n, cuDoubleComplex* A, int lda, cuFloatComplex  *SA , int * INFO , double RMAX , int ldsa )
{
    double mRMAX = - RMAX;
    int tx = threadIdx.x ;
    int ty = threadIdx.y ;
    int ind = blockIdx.x*  dgemv_bs + tx ;

    /*
      Reverse Computation ...
      - Left
      - Triangle
      - Up
    */
    A+= lda*(n-1) ;
    SA+= ldsa*(n-1) ;

    __shared__ cuDoubleComplex la   [dgemv_bs][dgemv_bs+1];

    A += ind;
    SA += ind;

    A-= ty * lda  ;
    SA-= ty * ldsa  ;

    cuDoubleComplex temp ;
    int break_d  = (n / dgemv_bs -   blockIdx.x-1 )* dgemv_bs ;

    for(int  i=0; i<break_d; i += dgemv_bs ){
#pragma unroll 8
        for(int j=0; j < dgemv_bs ; j+=4){
            //   la[tx][ty+j] = A[-j*lda] ;
            temp =  A[-j*lda] ;
            if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
                )
                {
                    flag = 1;
                }
            SA[-j*ldsa] = cuComplexDoubleToFloat(temp);
        }
        A-=lda* dgemv_bs ;
        SA-=ldsa* dgemv_bs ;
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
    for(int j =0; j<dgemv_bs; j+=4){
        temp =   la[tx][31-ty-j];
        if( (cuCreal(temp) < mRMAX) || (cuCreal(temp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
            || (cuCimag(temp) < mRMAX) || (cuCimag(temp) > RMAX) 
#endif
            )
            {
                flag = 1;
            }
        SA[- j*ldsa] = cuComplexDoubleToFloat(temp);
    }

    //la[tx][ty] = flag ;
    //
    //__syncthreads();
    //
    //if( ty == 0 ) {
    //    // INFO[0] = flag + la[tx][1] + la[tx][2] + la[tx][3] ;
    //}
}


extern "C" void
mzlat2c(char uplo, int m, cuDoubleComplex *A, int lda, cuFloatComplex *Y, int LDSA, int *INFO)
{
    /*
      Note:
      The UPLO = 'U' Version can be optimized more.
    */
    double RMAX = (double)lapackf77_slamch("O");
    int blocks;
    if (m % dgemv_bs==0)
        blocks = m/ dgemv_bs;
    else
        blocks = m/ dgemv_bs + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(32, 4, 1);

    if( m % dgemv_bs == 0 ) {
        if( uplo == 'L' || uplo == 'l'){
            l_zlat2c_special <<<grid, threads>>> (m, A, lda, Y ,INFO ,  RMAX , LDSA  );
        }
        else{
            u_zlat2c_special <<<grid, threads>>> (m, A, lda,  Y , INFO , RMAX , LDSA  );
        }

    }
    else{
        int  m_full_block = (m - m % 32 ) /32 ;
        int  m_mod_32 = m%32 ;
        if( uplo == 'L' || uplo == 'l'){
            l_zlat2c_generic <<<grid, threads>>> (m, A, lda, Y , m_full_block , m_mod_32 , INFO , RMAX , LDSA );
        }
        else{
            u_zlat2c_generic <<<grid, threads>>> (m, A, lda, Y , m_full_block , m_mod_32 , INFO , RMAX , LDSA );
        }
    }
}


/*
  Interface ..................................
  Reproduced from dlansy routines...
  How to deliver the info.
*/
extern "C" void
magmablas_zlat2c(char uplo, int n, cuDoubleComplex *A, int lda,
                 cuFloatComplex *SA, int LDSA, int *INFO)
{
    /*
      The routine converts a COMPLEX_16 triangular
      matrix A to COMPLEX triangular matrix SA.
    */
    *INFO = 0;
    mzlat2c( uplo, n, A, lda, SA, LDSA, INFO );
}

///////////////////////////////////////////////////////////////////////////////////////////
#else
///////////////////////////////////////////////////////////////////////////////////////////
/*------------------------------------------ UPLO = 'L' ----------------------------------*/

__global__ void
l_zlat2c_special (int n, cuDoubleComplex* A, int lda,  cuFloatComplex *SA , int *INFO, 
                  double RMAX, int ldsa )
{
    int tx  = threadIdx.x ;
    int ty  = threadIdx.y ;
    int ind = blockIdx.x*  dgemv_bs + tx ;

    A  += ind;
    SA += ind;
    A  += ty * lda;
    SA += ty * ldsa;

    int break_d  =   (blockIdx.x+1)* dgemv_bs ;
    for(int  i=0; i<break_d; i += dgemv_bs ){

        #pragma unroll 8
        for(int j=0; j < dgemv_bs ; j+=4)
            SA[j*ldsa] = cuComplexDoubleToFloat(A[j*lda]);

        A  += lda *dgemv_bs ;
        SA += ldsa*dgemv_bs ;
    }
}


__global__ void
l_zlat2c_generic(int n, cuDoubleComplex* A, int lda, cuFloatComplex *SA , int m_full_block,
                 int m_mod_32 , int *INFO , double RMAX , int ldsa)
{
    int tx = threadIdx.x ;
    int ty = threadIdx.y ;

    int ind = blockIdx.x*  dgemv_bs + tx ;

    if( blockIdx.x == m_full_block ) {
        /************************************************************************
   -- Last block --
   -- We will do something unusual here
   -- For sufficiently large matrix the overhead will be very low
        *************************************************************************/
        if  ( tx < m_mod_32 ){
            A+= ( blockIdx.x * dgemv_bs + tx ) ;
            SA+= ( blockIdx.x * dgemv_bs + tx ) ;
        }
        else{
            A+= ( blockIdx.x * dgemv_bs + m_mod_32 -1) ;
            SA+= ( blockIdx.x * dgemv_bs + m_mod_32 -1) ;
        }
        A+= ty * lda  ;
        SA+= ty * ldsa  ;
        int break_d  =   blockIdx.x* dgemv_bs ;

        /*----------------------------
          Go Right
          -------------------------------*/

        for(int  i=0; i<break_d; i += dgemv_bs ){
            #pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4)
                SA[j*ldsa] = cuComplexDoubleToFloat(A[j*lda]);
            
            A  += lda *dgemv_bs ;
            SA += ldsa*dgemv_bs ;
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
            for(j =0;j<=count;j++)
                SA[j*ldsa] = cuComplexDoubleToFloat(A[j*lda]);
            
            A  += (tx)*lda;
            SA += (tx)*ldsa;
            count = 1 ;
            for(;j<m_mod_32;j++){
                SA[count] = cuComplexDoubleToFloat(A[count]);
                count++;
            }
        }
        else{
        }

        __syncthreads();
    }

    else{
        /* **************************************
         -- All the blocks but the last one --
           ************************************** */
        A += ind;
        SA += ind;
        A+= ty * lda  ;
        SA+= ty * ldsa  ;
        int break_d  =   (blockIdx.x+1)* dgemv_bs ;

        /*----------------------------
          Go Right
          -------------------------------*/
        for(int  i=0; i<break_d; i += dgemv_bs ){
            #pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4)
                SA[j*ldsa] = cuComplexDoubleToFloat(A[j*lda]);
            
            A  += lda *dgemv_bs ;
            SA += ldsa*dgemv_bs ;
        }
   }
}

/*- ---------------------------------------------- UPLO = 'U' ----------------------------------*/
/* Generic Case*/

__global__ void
u_zlat2c_generic(int n, cuDoubleComplex* A, int lda, cuFloatComplex *SA, int m_full_block,
                 int m_mod_32 , int *INFO , double RMAX , int ldsa)
{
    int tx = threadIdx.x ;
    int ty = threadIdx.y ;

    int ind = blockIdx.x*  dgemv_bs + tx ;

    int blockIdxx =  blockIdx.x ;
    if( blockIdx.x == m_full_block ) {

        /************************************************************************
   -- Last block --
   -- We will do something unusual here
   -- For sufficiently large matrix the overhead will be very low
        *************************************************************************/
        ind =  tx ;
        A+= lda*(n-1) ;
        SA+= ldsa*(n-1) ;


        if  ( tx < m_mod_32 ){
            A+= (  tx ) ;
            SA+= (  tx ) ;
        }
        else{
            A+= (  m_mod_32 -1) ;
            SA+= (  m_mod_32 -1) ;
        }
        A-= ty * lda  ;
        SA-= ty * ldsa  ;
        int break_d  =   (blockIdx.x)* dgemv_bs ;

        /*----------------------------
          Go Right
          -------------------------------*/
        for(int  i=0; i<break_d; i += dgemv_bs ){
            #pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4)
                SA[-j*ldsa] = cuComplexDoubleToFloat(A[-j*lda]);
            
            A  -= lda *dgemv_bs ;
            SA -= ldsa*dgemv_bs ;
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
            for(j =0;j<count;j++)
                SA[-j*ldsa] = cuComplexDoubleToFloat(A[-j*lda]);
            
            A-=(count-1)*lda;
            SA-=(count-1)*ldsa;
            count = 1 ;
            for(;j<m_mod_32;j++){
                SA[-count] = cuComplexDoubleToFloat(A[-count]);
                count++;
            }
        }
        else{
        }
    }

    else{
        /* **************************************
           -- All the blocks but the last one --
           -- By the way this code can be optimized more.
           **************************************        */
        ind = blockIdx.x *  dgemv_bs + tx + m_mod_32 ;
        A+= lda*(n-1)  ;
        SA+= ldsa*(n-1)  ;

        A += ind;
        SA += ind;
        A-= ty * lda  ;
        SA-= ty * ldsa  ;

        int break_d  = (n / dgemv_bs -   blockIdxx )* dgemv_bs ;
        /*----------------------------
          Go Left
          -------------------------------*/
        for(int  i=0; i<break_d; i += dgemv_bs ){
            #pragma unroll 8
            for(int j=0; j < dgemv_bs ; j+=4)
                SA[-j*ldsa] = cuComplexDoubleToFloat(A[-j*lda]);
            
            A-=lda* dgemv_bs ;
            SA-=ldsa* dgemv_bs ;
        }
    }
}


/*- ---------------------------------------------- UPLO = 'U' ----------------------------------*/
/*Good Dimension*/
__global__ void
u_zlat2c_special (int n, cuDoubleComplex* A, int lda, cuFloatComplex  *SA , 
                  int * INFO , double RMAX , int ldsa )
{
    int tx = threadIdx.x ;
    int ty = threadIdx.y ;
    int ind = blockIdx.x*  dgemv_bs + tx ;

    /*
      Reverse Computation ...
      - Left
      - Triangle
      - Up
    */
    A+= lda*(n-1) ;
    SA+= ldsa*(n-1) ;

    A += ind;
    SA += ind;

    A-= ty * lda  ;
    SA-= ty * ldsa  ;

    int break_d  = (n / dgemv_bs -   blockIdx.x )* dgemv_bs ;

    for(int  i=0; i<break_d; i += dgemv_bs ){
        #pragma unroll 8
        for(int j=0; j < dgemv_bs ; j+=4)
            SA[-j*ldsa] = cuComplexDoubleToFloat(A[-j*lda]);
        
        A-=lda* dgemv_bs ;
        SA-=ldsa* dgemv_bs ;
    }
}


extern "C" void
mzlat2c(char uplo, int m, cuDoubleComplex *A, int lda, cuFloatComplex *Y, int LDSA, int *INFO)
{
    /*
      Note:
      The UPLO = 'U' Version can be optimized more.
    */
    double RMAX = (double)lapackf77_slamch("O");
    int blocks;
    if (m % dgemv_bs==0)
        blocks = m/ dgemv_bs;
    else
        blocks = m/ dgemv_bs + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(32, 4, 1);

    if( m % dgemv_bs == 0 ) {
        if( uplo == 'L' || uplo == 'l'){
            l_zlat2c_special <<<grid, threads>>> (m, A, lda, Y ,INFO ,  RMAX , LDSA  );
        }
        else{
            u_zlat2c_special <<<grid, threads>>> (m, A, lda,  Y , INFO , RMAX , LDSA  );
        }

    }
    else{
        int  m_full_block = (m - m % 32 ) /32 ;
        int  m_mod_32 = m%32 ;
        if( uplo == 'L' || uplo == 'l'){
            l_zlat2c_generic <<<grid, threads>>> (m, A, lda, Y , m_full_block , m_mod_32 , INFO , RMAX , LDSA );
        }
        else{
            u_zlat2c_generic <<<grid, threads>>> (m, A, lda, Y , m_full_block , m_mod_32 , INFO , RMAX , LDSA );
        }
    }
}


/*
  Interface ..................................
  Reproduced from dlansy routines...
  How to deliver the info.
*/
extern "C" void
magmablas_zlat2c(char uplo, int n, cuDoubleComplex *A, int lda,
                 cuFloatComplex *SA, int LDSA, int *INFO)
{
    /*
      The routine converts a COMPLEX_16 triangular
      matrix A (along with its block diagonal entries) to a COMPLEX 
      triangular matrix SA (along with its block diagonal entries).
    */
    *INFO = 0;
    mzlat2c( uplo, n, A, lda, SA, LDSA, INFO );
    /*
      int val = cublasIdamax(n,WORK,1);
      double retVal[1];
      cublasGetMatrix( 1, 1, sizeof( double ), WORK+val-1, 1, retVal, 1 ) ;
      return retVal[0];
    */
}

#endif /*  (!defined(PRECISION_z)) || (GPUSHMEM >= 200) */
