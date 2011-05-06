/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> z
       
*/
#include "common_magma.h"

#define magmablas_zgemv_tesla magmablas_zgemv

extern "C" void
magmablas_zgemv_tesla(char trans, int m, int n, 
                      cuDoubleComplex alpha, cuDoubleComplex *A, int lda, 
                                             cuDoubleComplex *x, int incx, 
                      cuDoubleComplex beta,  cuDoubleComplex *y, int incy) 
{
    cublasZgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
