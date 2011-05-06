/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated ds

*/
#include "common_magma.h"

static __global__ void 
slag2d_generic(int M, int N, 
               float *SA, int LDSA , 
               double *A , int LDA ) 
{ 
    int ibx = blockIdx.x * 64;

    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int idt = ty * 16 + tx;
	
    if( (ibx+idt) >= M ){
        SA += (M-1);
        A  += (M-1);
    }
    else{
        SA += ibx+idt;
        A  += ibx+idt;
    }
    const float * SAend = SA+LDSA*N;
    double Ap[1]={ (double)(SA[0]) };
    do {
        SA  += LDSA;
        A[0] = Ap[0];
        Ap[0]= (double)(SA[0]);
        A   += LDA;

    } while (SA < SAend);

    A[0] = Ap[0];
}

static __global__ void 
slag2d_special(int M, int N, 
               float *SA, int LDSA, 
               double *A , int LDA ) 
{ 
    int ibx = blockIdx.x * 64;

    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int idt = ty * 16 + tx;
	
    if( (ibx+idt) >= M ){
        SA += (M-1);
        A  += (M-1);
    }
    else{
        SA += ibx+idt;
        A  += ibx+idt;
    }
    double Ap[1] = { (double)(SA[0]) };
    A[0] = Ap[0];
}

static void 
magmablas_slag2d_64_64_16_4_v2( int M, int N, 
                                float *SA, int LDSA , 
                                double *A , int LDA )
{
    if( M == 0 || N==0 ) {
        printf("One of the dimension is ZERO\n");
        exit(-1);
    }
    dim3 threads( 16, 4 );
    dim3 grid(M/64+(M%64!=0),1);
    if( N > 1 ) {
        slag2d_generic<<< grid, threads >>> (  M, N, SA, LDSA, A, LDA ) ;
    }
    else{
        slag2d_special<<< grid, threads >>> (  M, N, SA, LDSA, A, LDA ) ;
    }
}

extern "C" void 
magmablas_slag2d(int M, int N, float *SA, int LDSA, double *A, int LDA, int *INFO)
{
/*
  Purpose
  =======

  CLAG2Z converts a SINGLE PRECISION matrix, SA, to a DOUBLE
  PRECISION matrix, A.

  Note that while it is possible to overflow while converting
  from double to single, it is not possible to overflow when
  converting from single to double.

  This is an auxiliary routine so there is no argument checking.

  Arguments
  =========

  M       (input) INTEGER
          The number of lines of the matrix A.  M >= 0.

  N       (input) INTEGER
          The number of columns of the matrix A.  N >= 0.

  SA      (input) REAL array, dimension (LDSA,N)
          On entry, the M-by-N coefficient matrix SA.

  LDSA    (input) INTEGER
          The leading dimension of the array SA.  LDSA >= max(1,M).

  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
          On exit, the M-by-N coefficient matrix A.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).

  INFO    (output) INTEGER
          = 0:  successful exit
  =========
*/
    *INFO = 0;
    magmablas_slag2d_64_64_16_4_v2( M, N, SA, LDSA, A, LDA ) ;
}	
