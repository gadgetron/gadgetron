/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> c

*/
#include "common_magma.h"
#define PRECISION_z

/*The version for tesla can be found in zhemv_tesla.cu */
#if (GPUSHMEM >= 200)

#define magmablas_zhemv_200 magmablas_zhemv

#define zhemv_bs         64
#define thread_x         64
#define thread_y          4
#define bank_shift       33
#define quarter_thread_x 16
#define half_thread_x    32

/*******************************************************************************
 *     Functions for each specific cases - Lower case
 */

__global__ void
magmablas_zhemv_200_L_special( magma_int_t n, cuDoubleComplex alpha,
                               cuDoubleComplex *A, magma_int_t lda,
                               cuDoubleComplex *x, magma_int_t incx,
                               cuDoubleComplex  beta,
                               cuDoubleComplex *y, magma_int_t incy,
                               cuDoubleComplex *WC)
{
    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;
    magma_int_t blkc = blockIdx.x ;

    cuDoubleComplex res  = MAGMA_Z_ZERO;
    cuDoubleComplex res_ = MAGMA_Z_ZERO;
    cuDoubleComplex res1 = MAGMA_Z_ZERO;

    __shared__ cuDoubleComplex la   [quarter_thread_x][thread_x+2];
    __shared__ cuDoubleComplex buff [thread_x];

    cuDoubleComplex tr[4];
    cuDoubleComplex b[4];

    magma_int_t break_d   =  thread_x * blkc;
    const magma_int_t td  = (thread_x * ty ) + tx;
    magma_int_t       tx_ = td % half_thread_x;
    magma_int_t       ty_ = td / half_thread_x;

    WC +=  break_d + tx;
    x  += (break_d + tx ) * incx;
    A  +=  break_d * (lda+1);
    A  += ty_* lda + tx_ ;

    if( ty == 0 ){
        buff[tx] = x[0];
    } // obtain the vector x store in buff;

    tx = tx_ ; ty = ty_ ;

    #pragma unroll
    for(magma_int_t j =0; j<half_thread_x; j +=8)
        la[0][ bank_shift * (ty_+j) + tx_] =  A[ j * lda];
    __syncthreads();

    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(ty_ * 4 + 4)  ; i++){
        if ( i < tx_ )   {
            la[0][bank_shift * tx_ + i] = cuConj( la[0][ i * bank_shift + tx_] ) ;
        }
        else
            la[0][bank_shift * tx_ + i] = la[0][ bank_shift * tx_ + i]  ;
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res+= cuConj( la[0][bank_shift * tx_ + j + ty_ * 4] ) * buff[j + ty_ * 4];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();

    if( ty_== 0 )
      res1 = la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
        +    la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
        +    la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
        +    la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    else
        {
            MAGMA_Z_SET2REAL(res1,0);
        }
    __syncthreads();


    MAGMA_Z_SET2REAL(res, 0) ;

    A+= half_thread_x + half_thread_x *lda ;

    #pragma unroll
    for(magma_int_t j =0; j<half_thread_x; j+=8)
        la[0][bank_shift*(ty_+j)+tx_] = A[ j * lda];
    __syncthreads();

    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(4+ty_*4) ; i++){
        if ( i < tx_ )   {
            la[0][bank_shift*tx_+i] = cuConj( la[0][bank_shift*i+tx_] ) ;
        }
        else
            la[0][bank_shift*tx_+i] = la[0][bank_shift*tx_+i]  ;
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res+= cuConj( la[0][bank_shift*tx_+j+ty_*4] ) * buff[half_thread_x + j + 4 * ty_];
    __syncthreads();
    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();

    cuDoubleComplex res2;
    MAGMA_Z_SET2REAL(res2,0);
    if( ty_== 1 )
        res2 = la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
          +    la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
          +    la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
          +    la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    else
    {
        MAGMA_Z_SET2REAL(res2,0);
    }
    __syncthreads();

    MAGMA_Z_SET2REAL(res,0);

    A-=half_thread_x *lda ;

    MAGMA_Z_SET2REAL(res_,0);

    #pragma unroll
    for(magma_int_t j=0; j<half_thread_x; j+=8)
        tr[j/8] = A[ j * lda];

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++){
        res += tr[j] * buff[ j*8 + ty_];
        la[0][bank_shift*(ty_+j*8)+tx_] = tr[j];
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res_+= cuConj(la[0][bank_shift*tx_+j+ty_*4]) * buff[half_thread_x +j+ty_*4];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();
    if( ty_ == 1 )
        res2 = res2 
            +  la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    else
        {
            MAGMA_Z_SET2REAL(res2,0);
        }
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res_ ;
    __syncthreads();
    if( ty_ == 0 ) {
        res1 = res1
            +  la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    }
    else
        {
            MAGMA_Z_SET2REAL(res1,0);
        }
    A-=half_thread_x;

    __syncthreads();
    tx = threadIdx.x ;
    ty = threadIdx.y ;

    if( ty_ == 0  && ty == 0  )
        res = res1 ;
    else if( ty_ == 1  && ty == 0  )
        res = res2 ;
    else
        {
            MAGMA_Z_SET2REAL(res,0);
        }

    A-=ty_* lda  ;
    A-=tx_;

    A= A - lda * blkc * thread_x;
    x= x - blkc * thread_x  *incx  ;

    x= x- tx*incx;

    A+=4 * ty* lda  ;
    A+=tx;

    magma_int_t wc_c = 0 ;
    magma_int_t count = 0 ;

    tx_ = td % quarter_thread_x ;
    ty_ = td / quarter_thread_x ;

    WC-=tx ;
    WC+=tx_;

    if( blkc * thread_x >=thread_x)
        #pragma unroll
        for(magma_int_t i=0; i<thread_x; i += thread_x )
        {
            MAGMA_Z_SET2REAL(res_,0);
            count++;

            #pragma unroll
            for( magma_int_t k=0;k<4;k++)
            {

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                    tr[j] = A[j*lda];

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                {
                    res += tr[j] * x[ quarter_thread_x * k + ty * 4 + j];
                    la[( j + ty * 4)][tx] = cuConj(tr[j]) * buff[tx];
                }
                __syncthreads();


                MAGMA_Z_SET2REAL(res_,0);

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                {
                    res_+=la[tx_][ty_*4+j] ;
                }
                b[k] = res_ ;
                __syncthreads();

                A += lda * quarter_thread_x ;
            }

            #pragma unroll
            for(magma_int_t k=0; k < 4 ; k++){
                la[tx_][ty_+quarter_thread_x*k]= b[k] ;
            }
            __syncthreads();
            if( ty_ < 4 ) {
                magma_int_t k = ty_*quarter_thread_x;
                res_ = la[tx_][0+k] + la[tx_][1+k]
                    +  la[tx_][2+k] + la[tx_][3+k]
                    +  la[tx_][4+k] + la[tx_][5+k]
                    +  la[tx_][6+k] + la[tx_][7+k]
                    +  la[tx_][8+k] + la[tx_][9+k]
                    +  la[tx_][10+k]+ la[tx_][11+k]
                    +  la[tx_][12+k]+ la[tx_][13+k]
                    +  la[tx_][14+k]+ la[tx_][15+k];
                WC[k + wc_c*lda ] =   res_;
            }

            wc_c++;
            __syncthreads();

        }

    for(magma_int_t  i=thread_x; i< (blkc * thread_x); i += thread_x )
    {
        MAGMA_Z_SET2REAL(res_,0);
        count++;

        #pragma unroll
        for( magma_int_t k=0;k<4;k++)
	{
            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
                tr[j] = A[j*lda] ;

            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
            {
                res += tr[j] * x[ i + quarter_thread_x*k + ty*4+(j)];
                la[( j + ty * 4)][tx] = cuConj( tr[j] )* buff[tx];
            }
            __syncthreads();

            MAGMA_Z_SET2REAL(res_,0);

            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
                res_+=la[tx_][ty_*4+j] ;

            b[k] = res_ ;
            __syncthreads();

            A += lda * quarter_thread_x ;
        }

        #pragma unroll
        for(magma_int_t k=0; k < 4 ; k++){
            la[tx_][ty_+quarter_thread_x*k]= b[k] ;
        }
        __syncthreads();
        if( ty_ < 4 ) {
            magma_int_t k = ty_*quarter_thread_x;
            res_ = la[tx_][0+k] + la[tx_][1+k]
                +  la[tx_][2+k] + la[tx_][3+k]
                +  la[tx_][4+k] + la[tx_][5+k]
                +  la[tx_][6+k] + la[tx_][7+k]
                +  la[tx_][8+k] + la[tx_][9+k]
                +  la[tx_][10+k]+ la[tx_][11+k]
                +  la[tx_][12+k]+ la[tx_][13+k]
                +  la[tx_][14+k]+ la[tx_][15+k];
            WC[k + wc_c*lda ] =   res_;
        }

        wc_c++;
        __syncthreads();
    }

    WC+=tx ;
    WC-=tx_;

    la[ty][tx]= res ;
    __syncthreads();
    if( ty == 0 ) {
        res = la[0][tx]+ la[1][tx]
            + la[2][tx]+ la[3][tx];
        WC[0+lda*(blkc)  ] =  res;
    }
}

/**************************************************************
 *    Lower case for generic sizes
 */
__global__ void
magmablas_zhemv_200_L_generic(magma_int_t n, cuDoubleComplex alpha,
                              cuDoubleComplex *A, magma_int_t lda,
                              cuDoubleComplex *x, magma_int_t incx,
                              cuDoubleComplex beta,
                              cuDoubleComplex *y, magma_int_t incy,
                              cuDoubleComplex *WC,
                              magma_int_t m_mod_thread_x)
{
    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;
    magma_int_t blkc = blockIdx.x ;

    cuDoubleComplex res  = MAGMA_Z_ZERO;
    cuDoubleComplex res_ = MAGMA_Z_ZERO;
    cuDoubleComplex res1 = MAGMA_Z_ZERO;

    __shared__ cuDoubleComplex la   [quarter_thread_x][thread_x+2];
    __shared__ cuDoubleComplex buff [thread_x];
    __shared__ cuDoubleComplex buff2[thread_x];

    cuDoubleComplex tr[4];
    cuDoubleComplex b[8];

    magma_int_t break_d   =  thread_x * blkc;
    const magma_int_t td  = (thread_x * ty ) + tx;
    magma_int_t       tx_ = td % half_thread_x;
    magma_int_t       ty_ = td / half_thread_x;

    WC+=  break_d + tx;
    x += (break_d + tx ) * incx;
    A +=  break_d * (lda+1);
    A += lda * ty_;

    magma_int_t trackA ;
    if( blkc == ( gridDim.x - 1 ) ) {
        if( ty == 0 ){
            if( tx > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(buff[tx],0);
            }
            else
                buff[tx]  = x[0];
        }
        if ( tx_ > m_mod_thread_x )
            trackA=m_mod_thread_x;
        else
            trackA=tx_;
        A += trackA ;
    }
    else {
        if( ty == 0 ){
            buff[tx]  = x[0];
        }
        trackA = tx_;
        A += trackA ;
    }

    // Somehow merging these two if - else creates problem
    // It could be a potential bug -- from synchronization or from cuda or compiler
    if( blkc == ( gridDim.x - 1 ) ) {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            if( ( ty_ + j ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(la[0][bank_shift*(ty_+j)+tx_], 9999);
            }
            else
                la[0][bank_shift*(ty_+j)+tx_] =  A[ j * lda];
        }
        A-=trackA;
    }
    else {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            la[0][bank_shift*(ty_+j)+tx_] = A[ j * lda];
        }
    }
    tx = tx_ ;
    ty = ty_ ;
    __syncthreads();

    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(ty_*4+4)  ; i++){
        if ( i < tx_ )   {
            la[0][bank_shift*tx_+i] = cuConj(la[0][i*bank_shift+tx_]) ;
        }
        else
            la[0][bank_shift*tx_+i] = la[0][bank_shift*tx_+i]  ;
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res += cuConj(la[0][bank_shift*tx_+j+ty_*4])* buff[j+ty_*4];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();
    if( ty_== 0 )
        res1 = la[0][tx_*bank_shift+0] 
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    else
    {
        MAGMA_Z_SET2REAL(res1,0);
    }
    __syncthreads();


    MAGMA_Z_SET2REAL(res,0);

    if( blkc == ( gridDim.x - 1 ) ) {
        if ( (tx_+half_thread_x) > m_mod_thread_x )
            trackA = m_mod_thread_x;
        else
            trackA = tx_ + half_thread_x;
        A+= trackA+half_thread_x*lda ;

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            if( ( ty_ + j+half_thread_x ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(la[0][bank_shift*(ty_+j)+tx_], 99999);
            }
            else
                la[0][bank_shift*(ty_+j)+tx_] =  A[ j * lda];
        }

        A-= trackA+half_thread_x*lda ;
        A+=tx_ ;
        A+= half_thread_x + half_thread_x *lda ;
    }
    else {
        A+= half_thread_x + half_thread_x *lda ;

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            la[0][bank_shift*(ty_+j)+tx_] = A[ j * lda];
        }
    }

    __syncthreads();
    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(4+ty_*4) ; i++){
        if ( i < tx_ )   {
            la[0][bank_shift*tx_+i] = cuConj(la[0][bank_shift*i+tx_]) ;
        }
        else
            la[0][bank_shift*tx_+i] = la[0][bank_shift*tx_+i]  ;
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res+= cuConj(la[0][bank_shift*tx_+j+ty_*4]) * buff[half_thread_x + j + 4 * ty_];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();

    cuDoubleComplex res2;
    MAGMA_Z_SET2REAL(res2,0);
    if( ty_== 1 )
        res2 = la[0][tx_*bank_shift+0]
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    else
    {
        MAGMA_Z_SET2REAL(res2,0);
    }
    __syncthreads();

    MAGMA_Z_SET2REAL(res,0);
    MAGMA_Z_SET2REAL(res_,0);

    A-=half_thread_x *lda ;
    if( blkc == ( gridDim.x - 1 ) ) {
        A-=tx_;
        if ( tx_ > m_mod_thread_x )
            trackA=m_mod_thread_x;
        else
            trackA=tx_;
        A+= trackA ;

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8)
            if( ( ty_ + j ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(tr[j/8], 99999);
            }
            else
                tr[j/8] = A[ j * lda];
        A-=trackA;
        A+=tx_;
    }
    else {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8)
            tr[j/8] = A[ j * lda];
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++){
        res+= tr[j] * buff[ j*8 + ty_];
        la[0][bank_shift*(ty_+j*8)+tx_] = tr[j];
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res_+= cuConj(la[0][bank_shift*tx_+j+ty_*4]) * buff[half_thread_x +j+ty_*4];
    __syncthreads();


    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();
    if( ty_ == 1 )
        res2 = res2
            +  la[0][tx_*bank_shift+0]
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    else
    {
        MAGMA_Z_SET2REAL(res2,0);
    }
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res_ ;
    __syncthreads();

    if( ty_ == 0 ) {
        res1 = res1
            +  la[0][tx_*bank_shift+0]
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    }
    else
    {
        MAGMA_Z_SET2REAL(res1,0);
    }
    A-=half_thread_x;

    __syncthreads();
    tx = threadIdx.x ;
    ty = threadIdx.y ;

    if( ty_ == 0  && ty == 0  )
        res = res1 ;
    else if( ty_ == 1  && ty == 0  )
        res = res2 ;
    else
    {
        MAGMA_Z_SET2REAL(res,0);
    }

    A-=ty_* lda  ;
    A-=tx_;

    A= A - lda*break_d;
    x= x - break_d *incx  ;

    A+=4 * ty* lda  ;

    if( blkc  == ( gridDim.x - 1 ) ) {
        if(tx <= m_mod_thread_x )
            A+=tx;
        else
            A+=m_mod_thread_x;
    }
    else{
        A+=tx;
    }

    magma_int_t wc_c = 0 ;
    magma_int_t count = 0 ;

    tx_ = td % quarter_thread_x ;
    ty_ = td / quarter_thread_x ;

    WC-=tx ;
    WC+=tx_;

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        b[j] =  buff[ty_*4+j];

    if( break_d > 0)
        #pragma unroll
        for(magma_int_t  i=0; i< thread_x; i += thread_x ){
            MAGMA_Z_SET2REAL(res_,0);
            count++;
            if( ty== 0 ) {
                buff2[tx]  = x[i*incx];
            }
            __syncthreads();

            #pragma unroll
            for( magma_int_t k=0;k<4;k++){
                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                    tr[j] = A[j*lda] ;

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++){
                    res+=tr[j]*buff2[quarter_thread_x*k + ty*4+(j)];
                    la[( (j)+ty*4)][tx] = cuConj(tr[j]);
                }
                __syncthreads();

                MAGMA_Z_SET2REAL(res_, 0) ;

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                    res_+=la[tx_][ty_*4+j]* b[j] ;
                b[4+k] = res_ ;
                __syncthreads();
                A+=lda* quarter_thread_x ;
            }

            #pragma unroll
            for(magma_int_t k=0; k < 4 ; k++){
                la[tx_][ty_+quarter_thread_x*k]= b[4+k] ;
            }
            __syncthreads();

            if( ty_ < 4 ) {
                magma_int_t k = ty_*quarter_thread_x;
                res_ = la[tx_][0+k] + la[tx_][1+k] 
                    +  la[tx_][2+k] + la[tx_][3+k]
                    +  la[tx_][4+k] + la[tx_][5+k]
                    +  la[tx_][6+k] + la[tx_][7+k]
                    +  la[tx_][8+k] + la[tx_][9+k]
                    +  la[tx_][10+k]+ la[tx_][11+k]
                    +  la[tx_][12+k]+ la[tx_][13+k]
                    +  la[tx_][14+k]+ la[tx_][15+k];
                WC[k + wc_c*lda ] =   res_;
            }
            wc_c++;
            __syncthreads();
        }

    for(magma_int_t  i=thread_x; i<break_d; i += thread_x ){
        MAGMA_Z_SET2REAL(res_, 0) ;
        count++;
        if(ty == 0 )
            buff2[tx]  = x[i*incx];
        __syncthreads();

        #pragma unroll
        for( magma_int_t k=0;k<4;k++){
	    #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
                tr[j] = A[j*lda] ;
            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++){
                res+=tr[j]*buff2[quarter_thread_x*k + ty*4+(j)];
                la[( (j)+ty*4)][tx] = cuConj(tr[j]);
	    }
	    __syncthreads();

            MAGMA_Z_SET2REAL(res_, 0) ;

	    #pragma unroll
	    for(magma_int_t j=0; j < 4 ; j++)
                res_+=la[tx_][ty_*4+j]* b[j] ;
	    b[4+k] = res_ ;
	    __syncthreads();
	    A+=lda* quarter_thread_x ;
        }

        #pragma unroll
        for(magma_int_t k=0; k < 4 ; k++){
            la[tx_][ty_+quarter_thread_x*k]= b[4+k] ;
        }
        __syncthreads();

        if( ty_ < 4 ) {
            magma_int_t k = ty_*quarter_thread_x;
            res_ = la[tx_][0+k] + la[tx_][1+k] 
                +  la[tx_][2+k] + la[tx_][3+k]
                +  la[tx_][4+k] + la[tx_][5+k]
                +  la[tx_][6+k] + la[tx_][7+k]
                +  la[tx_][8+k] + la[tx_][9+k]
                +  la[tx_][10+k]+ la[tx_][11+k]
                +  la[tx_][12+k]+ la[tx_][13+k]
                +  la[tx_][14+k]+ la[tx_][15+k];
            WC[k + wc_c*lda ] =   res_;
        }
        wc_c++;
        __syncthreads();
    }


    WC+=tx ;
    WC-=tx_;
    la[ty][tx]= res ;
    __syncthreads();

    if( ty == 0 ) {
        res=la[0][tx]+ la[1][tx]+ la[2][tx]+ la[3][tx] ;
        WC[0+lda*(blkc)] = res;
    }
}

__global__ void
magmablas_zhemv_200_L_update(magma_int_t n, cuDoubleComplex alpha,
                         cuDoubleComplex* A, magma_int_t lda,
                         cuDoubleComplex *x, magma_int_t incx,
                         cuDoubleComplex beta,
                         cuDoubleComplex *y, magma_int_t incy,
                         cuDoubleComplex *WC )
{
    magma_int_t i;
    magma_int_t tx  = threadIdx.x ;
    magma_int_t ind = blockIdx.x * thread_x + tx ;
    cuDoubleComplex Ca;

    MAGMA_Z_SET2REAL(Ca, 0) ;
    WC+= ind + lda * blockIdx.x;

    for(i = blockIdx.x*thread_x; i<n; i+=thread_x){
        Ca += WC[0] ;
        WC += thread_x;
    }
    if( ind < n )
        y[ind * incy] = beta * y[ind * incy]  + alpha * Ca ;
}


extern "C"
void magmablas_zhemv_200_L(magma_int_t m, cuDoubleComplex alpha,
                           cuDoubleComplex *A, magma_int_t lda,
                           cuDoubleComplex *X, magma_int_t incx,
                           cuDoubleComplex beta,
                           cuDoubleComplex *Y, magma_int_t incy,
                           cuDoubleComplex *dC_work)
{
    magma_int_t blocks;

    if (m % zhemv_bs==0)
        blocks = m / zhemv_bs;
    else
        blocks = m / zhemv_bs + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(thread_x, thread_y, 1);
    dim3 threads_u(zhemv_bs, 1, 1);

    /*
     * If matrix size is multiple of zhemv_bs, we use a specific code.
     * otherwise, we call the generic case.
     */
    if(m % zhemv_bs == 0 ) {
        magmablas_zhemv_200_L_special <<<grid, threads>>>(
            m, alpha, A, lda, X, incx, beta, Y, incy, dC_work);
    }
    else{
        magma_int_t m_mod_thread_x = m%zhemv_bs - 1;
        magmablas_zhemv_200_L_generic <<<grid, threads>>> (
            m, alpha, A, lda, X, incx ,beta, Y, incy, dC_work, m_mod_thread_x);
    }

    magmablas_zhemv_200_L_update<<<grid, threads_u>>>(
        m, alpha, A, lda, X, incx, beta, Y, incy, dC_work);
}

/*************************************************************************

    Purpose
    =======

    magmablas_zhemv  performs the matrix-vector operation on fermi:

       y := alpha*A*x + beta*y,

    where alpha and beta are scalars, x and y are n element vectors and
    A is an n by n hermitian matrix.

    Arguments
    ==========

    UPLO   - CHARACTER*1.
             On entry, UPLO specifies whether the upper or lower
             triangular part of the array A is to be referenced as
             follows:

                UPLO = 'U' or 'u'   Only the upper triangular part of A
                                    is to be referenced.

                UPLO = 'L' or 'l'   Only the lower triangular part of A
                                    is to be referenced.

             Unchanged on exit.

    N      - INTEGER.
             On entry, N specifies the order of the matrix A.
             N must be at least zero.
             Unchanged on exit.

    ALPHA  - COMPLEX*16      .
             On entry, ALPHA specifies the scalar alpha.
             Unchanged on exit.

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
             Before entry with  UPLO = 'U' or 'u', the leading n by n
             upper triangular part of the array A must contain the upper
             triangular part of the hermitian matrix and the strictly
             lower triangular part of A is not referenced.
             Before entry with UPLO = 'L' or 'l', the leading n by n
             lower triangular part of the array A must contain the lower
             triangular part of the hermitian matrix and the strictly
             upper triangular part of A is not referenced.
             Note that the imaginary parts of the diagonal elements need
             not be set and are assumed to be zero.
             Unchanged on exit.

    LDA    - INTEGER.
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program. LDA must be at least
             max( 1, n ).
             Unchanged on exit.
             It is recommended that lda is multiple of 16. Otherwise
             performance would be deteriorated as the memory accesses
             would not be fully coalescent.

    X      - COMPLEX*16       array of dimension at least
             ( 1 + ( n - 1 )*abs( INCX ) ).
             Before entry, the incremented array X must contain the n
             element vector x.
             Unchanged on exit.

    INCX   - INTEGER.
             On entry, INCX specifies the increment for the elements of
             X. INCX must not be zero.
             Unchanged on exit.

    BETA   - COMPLEX*16      .
             On entry, BETA specifies the scalar beta. When BETA is
             supplied as zero then Y need not be set on input.
             Unchanged on exit.

    Y      - COMPLEX*16       array of dimension at least
             ( 1 + ( n - 1 )*abs( INCY ) ).
             Before entry, the incremented array Y must contain the n
             element vector y. On exit, Y is overwritten by the updated
             vector y.

    INCY   - INTEGER.
             On entry, INCY specifies the increment for the elements of
             Y. INCY must not be zero.
             Unchanged on exit.

*/

extern "C"
magma_int_t
magmablas_zhemv_200( char uplo, magma_int_t n,
                     cuDoubleComplex alpha, 
                     cuDoubleComplex *A, magma_int_t lda,
                     cuDoubleComplex *X, magma_int_t incx,
                     cuDoubleComplex beta,  
                     cuDoubleComplex *Y, magma_int_t incy)
{
    char      uplo_[2] = {uplo, 0};
    long int  upper    = lapackf77_lsame(uplo_, "U");

    /*
     * Test the input parameters.
     */
    if ((! upper) && (! lapackf77_lsame(uplo_, "L"))) {
        return -1;
    } else if ( n < 0 ) {
        return -2;
    } else if ( lda < max(1,n) ) {
        return -5;
    } else if ( incx == 0 ) {
        return -7;
    } else if ( incy == 0 ) {
        return -10;
    }

    /*
     * Quick return if possible.
     */
    if ( (n == 0) || ( MAGMA_Z_EQUAL(alpha, MAGMA_Z_ZERO) && MAGMA_Z_EQUAL(beta, MAGMA_Z_ONE) ) )
        return MAGMA_SUCCESS;

    /* TODO: Upper case is not implemented in MAGMA */
    if ( upper )
        cublasZhemv(uplo, n, alpha, A, lda, X, incx, beta, Y, incy);
    else
    {
	cuDoubleComplex *dC_work;
	magma_int_t blocks    = n / thread_x + (n % thread_x != 0);
	magma_int_t workspace = lda * (blocks + 1);

        /* TODO: need to add a MAGMA context to handle workspaces */
	cublasAlloc( workspace, sizeof(cuDoubleComplex), (void**)&dC_work ) ;
        cublasGetError( ) ;

	magmablas_zhemv_200_L(n, alpha, A, lda, X, incx, beta, Y, incy, dC_work);

	cublasFree(dC_work);
        cublasGetError( ) ;
    }
    return MAGMA_SUCCESS;
}

#endif /* (GPUSHMEM >= 200) */
