/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated ds

*/
#include "common_magma.h"

#define PRECISION_d
#define blksize 64

__device__ int flag = 0; 

static __global__ void 
magmaint_dlag2s( magma_int_t M, magma_int_t N, 
                  const double *A, int lda, 
                  float *SA,       int ldsa, 
                  double RMAX ) 
{
    const double *Aend = A + lda*N;
    double tmp;
    double mRMAX = - RMAX;
    int    mym   = blockIdx.x * blksize + threadIdx.x;

    if ( mym < M ){
        A += mym;
        SA+= mym; 
        
        tmp = *A;
        for ( ; A < Aend; )
        {
            A  += lda;
            if( ((tmp) < mRMAX) || ((tmp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                || ((tmp) < mRMAX) || ((tmp) > RMAX) 
#endif
                )
            {
                flag = 1; 
            }
            *SA = (float)( tmp );
            tmp = *A;
            SA += ldsa;
        }
    }
}


extern "C" void 
magmablas_dlag2s( int M, int N , 
                  const double *A, int lda, 
                  float *SA,       int ldsa, 
                  magma_int_t *info ) 
{    
/*
  Note
  ====
	- We have to provide INFO at the end that dlag2s isn't doable now. 
	- Transfer a single value TO/FROM CPU/GPU
	- SLAMCH that's needed is called from underlying BLAS
	- Only used in iterative refinement
	- Do we want to provide this in the release?

  Purpose
  =======
  DLAG2S converts a DOUBLE PRECISION matrix A to a SINGLE PRECISION
  matrix SA.

  RMAX is the overflow for the SINGLE PRECISION arithmetic.
  DLAG2S checks that all the entries of A are between -RMAX and
  RMAX. If not the convertion is aborted and a flag is raised.

  This is an auxiliary routine so there is no argument checking.

  Arguments
  =========
  M       (input) INTEGER
          The number of lines of the matrix A.  M >= 0.

  N       (input) INTEGER
          The number of columns of the matrix A.  N >= 0.

  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the M-by-N coefficient matrix A.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).

  SA      (output) SINGLE PRECISION array, dimension (LDSA,N)
          On exit, if INFO=0, the M-by-N coefficient matrix SA; if
          INFO>0, the content of SA is unspecified.

  LDSA    (input) INTEGER
          The leading dimension of the array SA.  LDSA >= max(1,M).

  INFO    (output) INTEGER
          = 0:  successful exit.
          = 1:  an entry of the matrix A is greater than the SINGLE PRECISION
                overflow threshold, in this case, the content
                of SA in exit is unspecified.
  ===========================================================================  */

    double RMAX = (double)lapackf77_slamch("O");

    dim3 threads( blksize, 1, 1 );
    dim3 grid( (M+blksize-1)/blksize, 1, 1);
    flag = 0;
    magmaint_dlag2s<<< grid, threads >>>( M, N, A, lda, SA, ldsa, RMAX ) ; 
    *info = flag;
}
