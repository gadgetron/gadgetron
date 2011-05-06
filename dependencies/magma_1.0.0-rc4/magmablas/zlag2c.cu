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
#define blksize 64

__device__ int flag = 0; 

static __global__ void 
magmaint_zlag2c( magma_int_t M, magma_int_t N, 
                  const cuDoubleComplex *A, int lda, 
                  cuFloatComplex *SA,       int ldsa, 
                  double RMAX ) 
{
    const cuDoubleComplex *Aend = A + lda*N;
    cuDoubleComplex tmp;
    double mRMAX = - RMAX;
    int    mym   = blockIdx.x * blksize + threadIdx.x;

    if ( mym < M ){
        A += mym;
        SA+= mym; 
        
        tmp = *A;
        for ( ; A < Aend; )
        {
            A  += lda;
            if( (cuCreal(tmp) < mRMAX) || (cuCreal(tmp) > RMAX)
#if defined(PRECISION_z) || defined(PRECISION_c)
                || (cuCimag(tmp) < mRMAX) || (cuCimag(tmp) > RMAX) 
#endif
                )
            {
                flag = 1; 
            }
            *SA = cuComplexDoubleToFloat( tmp );
            tmp = *A;
            SA += ldsa;
        }
    }
}


extern "C" void 
magmablas_zlag2c( int M, int N , 
                  const cuDoubleComplex *A, int lda, 
                  cuFloatComplex *SA,       int ldsa, 
                  magma_int_t *info ) 
{    
/*
  Note
  ====
	- We have to provide INFO at the end that zlag2c isn't doable now. 
	- Transfer a single value TO/FROM CPU/GPU
	- SLAMCH that's needed is called from underlying BLAS
	- Only used in iterative refinement
	- Do we want to provide this in the release?

  Purpose
  =======
  ZLAG2C converts a COMPLEX_16 matrix A to a COMPLEX
  matrix SA.

  RMAX is the overflow for the COMPLEX arithmetic.
  ZLAG2C checks that all the entries of A are between -RMAX and
  RMAX. If not the convertion is aborted and a flag is raised.

  This is an auxiliary routine so there is no argument checking.

  Arguments
  =========
  M       (input) INTEGER
          The number of lines of the matrix A.  M >= 0.

  N       (input) INTEGER
          The number of columns of the matrix A.  N >= 0.

  A       (input) COMPLEX_16 array, dimension (LDA,N)
          On entry, the M-by-N coefficient matrix A.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).

  SA      (output) COMPLEX array, dimension (LDSA,N)
          On exit, if INFO=0, the M-by-N coefficient matrix SA; if
          INFO>0, the content of SA is unspecified.

  LDSA    (input) INTEGER
          The leading dimension of the array SA.  LDSA >= max(1,M).

  INFO    (output) INTEGER
          = 0:  successful exit.
          = 1:  an entry of the matrix A is greater than the COMPLEX
                overflow threshold, in this case, the content
                of SA in exit is unspecified.
  ===========================================================================  */

    double RMAX = (double)lapackf77_slamch("O");

    dim3 threads( blksize, 1, 1 );
    dim3 grid( (M+blksize-1)/blksize, 1, 1);
    flag = 0;
    magmaint_zlag2c<<< grid, threads >>>( M, N, A, lda, SA, ldsa, RMAX ) ; 
    *info = flag;
}
