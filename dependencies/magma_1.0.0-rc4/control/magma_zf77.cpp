/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include "magma.h"

/* 
 * typedef comming from fortran.h file provided in $CUDADIR/src directory
 * it will probably change with future release of cublas when they will use 64bits address
 */
typedef size_t devptr_t;

#undef REAL
#define COMPLEX

#ifdef ADD_

    #define MAGMA_ZGEBRD magma_zgebrd_
    #define MAGMA_ZGEHRD2 magma_zgehrd2_
    #define MAGMA_ZGELQF magma_zgelqf_
    #define MAGMA_ZGEQLF magma_zgeqlf_
    #define MAGMA_ZGEQRF magma_zgeqrf_
    #define MAGMA_ZGETRF magma_zgetrf_
    #define MAGMA_ZLAHR2 magma_zlahr2_
    #define MAGMA_ZLAHRU magma_zlahru_
    #define MAGMA_ZPOTRF magma_zpotrf_
    #define MAGMA_ZHETRD magma_zhetrd_
    
    #define MAGMA_ZUNMQR_GPU magma_zunmqr_gpu_
    #define MAGMA_ZGEQRF_GPU  magma_zgeqrf_gpu_
    #define MAGMA_ZGEQRF2_GPU magma_zgeqrf2_gpu_
    #define MAGMA_ZGEQRS_GPU magma_zgeqrs_gpu_
    #define MAGMA_ZGETRF_GPU magma_zgetrf_gpu_
    #define MAGMA_ZGETRS_GPU magma_zgetrs_gpu_
    #define MAGMA_ZGESV_GPU  magma_zgesv_gpu_
    #define MAGMA_ZLABRD_GPU magma_zlabrd_gpu_
    #define MAGMA_ZLARFB_GPU magma_zlarfb_gpu_
    #define MAGMA_ZPOTRF_GPU magma_zpotrf_gpu_
    #define MAGMA_ZPOTRS_GPU magma_zpotrs_gpu_

#elif defined (NOCHANGE)

    #define MAGMA_ZGEBRD magma_zgebrd
    #define MAGMA_ZGEHRD2 magma_zgehrd2
    #define MAGMA_ZGELQF magma_zgelqf
    #define MAGMA_ZGEQLF magma_zgeqlf
    #define MAGMA_ZGEQRF magma_zgeqrf
    #define MAGMA_ZGETRF magma_zgetrf
    #define MAGMA_ZLAHR2 magma_zlahr2
    #define MAGMA_ZLAHRU magma_zlahru
    #define MAGMA_ZPOTRF magma_zpotrf
    #define MAGMA_ZHETRD magma_zhetrd
    
    #define MAGMA_ZUNMQR_GPU magma_zunmqr_gpu
    #define MAGMA_ZGEQRF_GPU  magma_zgeqrf_gpu
    #define MAGMA_ZGEQRF2_GPU magma_zgeqrf2_gpu
    #define MAGMA_ZGEQRS_GPU magma_zgeqrs_gpu
    #define MAGMA_ZGETRF_GPU magma_zgetrf_gpu
    #define MAGMA_ZGETRS_GPU magma_zgetrs_gpu
    #define MAGMA_ZGESV_GPU  magma_zgesv_gpu
    #define MAGMA_ZLABRD_GPU magma_zlabrd_gpu
    #define MAGMA_ZLARFB_GPU magma_zlarfb_gpu
    #define MAGMA_ZPOTRF_GPU magma_zpotrf_gpu
    #define MAGMA_ZPOTRS_GPU magma_zpotrs_gpu

#endif

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
void MAGMA_ZGEBRD( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double *d, double *e, double2 *tauq,  double2 *taup, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_zgebrd( *m, *n, A, *lda, d, e, tauq,  taup, work, *lwork, info); 
}

void MAGMA_ZGEHRD2( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_zgehrd2( *n, *ilo, *ihi, A, *lda, tau, work, lwork, info); 
}

void MAGMA_ZGELQF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_zgelqf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_ZGEQLF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_zgeqlf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_ZGEQRF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_zgeqrf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_ZGETRF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    magma_zgetrf( *m, *n, A, *lda, ipiv, info); 
}

void MAGMA_ZLABRD_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                       double2 *a, magma_int_t *lda, devptr_t *da, magma_int_t *ldda, 
                       double *d, double *e, double2 *tauq, double2 *taup, 
                       double2 *x, magma_int_t *ldx, devptr_t *dx, magma_int_t *lddx,
                       double2 *y, magma_int_t *ldy, devptr_t *dy, magma_int_t *lddy)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*da);
    cuDoubleComplex *d_x = (cuDoubleComplex *)(*dx);
    cuDoubleComplex *d_y = (cuDoubleComplex *)(*dy);
    magma_zlabrd_gpu( *m, *n, *nb, 
                      a, *lda, d_a, *ldda, 
                      d, e, tauq, taup, 
                      x, *ldx, d_x, *lddx,
                      y, *ldy, d_y, *lddy); 
}

void MAGMA_ZLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, devptr_t *da, devptr_t *dv, double2 *a, magma_int_t *lda, double2 *tau, double2 *t, magma_int_t *ldt, double2 *y, magma_int_t *ldy)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*da);
    cuDoubleComplex *d_v = (cuDoubleComplex *)(*dv);
    magma_zlahr2( *m, *n, *nb, d_a, d_v, a, *lda, tau, t, *ldt, y, *ldy); 
}

void MAGMA_ZLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, double2 *a, magma_int_t *lda, devptr_t *da, double2 *y, double2 *v, double2 *t, devptr_t *dwork)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*da);
    cuDoubleComplex *d_w = (cuDoubleComplex *)(*dwork);
    magma_zlahru( *m, *n, *nb, a, *lda, d_a, y, v, t, d_w); 
}

void MAGMA_ZPOTRF( char *uplo, magma_int_t *n, double2 *A, magma_int_t *lda, magma_int_t *info)
{ 
    magma_zpotrf( uplo[0], *n, A, *lda, info); 
}

void MAGMA_ZHETRD( char *uplo, magma_int_t *n, double2 *A, magma_int_t *lda, double *d, double *e, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_zhetrd( uplo[0], *n, A, *lda, d, e, tau, work, *lwork, info); 
}


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMA_ZUNMQR_GPU(char *side, char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *a, magma_int_t *lda, double2 *tau, devptr_t *c, magma_int_t *ldc, double2 *work, magma_int_t *lwork, devptr_t *td, magma_int_t *nb, magma_int_t *info)
{ 
    cuDoubleComplex *d_a  = (cuDoubleComplex *)(*a);
    cuDoubleComplex *d_c  = (cuDoubleComplex *)(*c);
    cuDoubleComplex *d_td = (cuDoubleComplex *)(*td);
    magma_zunmqr_gpu(side[0], trans[0], *m, *n, *k, d_a, *lda, tau, d_c, *ldc, work, *lwork, d_td, *nb, info); 
}

void MAGMA_ZGEQRF2_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, double2 *tau, magma_int_t *info)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    magma_zgeqrf2_gpu( *m, *n, d_a, *lda, tau, info); 
}

void MAGMA_ZGEQRF_GPU(magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, double2 *tau, devptr_t *dwork, magma_int_t *info)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    cuDoubleComplex *d_w = (cuDoubleComplex *)(*dwork);
    magma_zgeqrf_gpu(*m, *n, d_a, *lda, tau, d_w, info); 
}

void MAGMA_ZGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
		       devptr_t *A, magma_int_t *lda, double2 *tau, devptr_t *td, 
		       devptr_t *c, magma_int_t *ldc, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    cuDoubleComplex *d_a  = (cuDoubleComplex *)(*A);
    cuDoubleComplex *d_c  = (cuDoubleComplex *)(*c);
    cuDoubleComplex *d_td = (cuDoubleComplex *)(*td);
    magma_zgeqrs_gpu( *m, *n, *nrhs, d_a, *lda, tau, d_td, d_c, *ldc, work, *lwork, info); 
}

void MAGMA_ZGETRF_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    magma_zgetrf_gpu( *m, *n, d_a, *lda, ipiv, info); 
}

void MAGMA_ZGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
  cuDoubleComplex *d_a = (cuDoubleComplex *)(uintptr_t)(*A);
  cuDoubleComplex *d_b = (cuDoubleComplex *)(uintptr_t)(*b);
  magma_zgetrs_gpu( trans[0], *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info); 
}

void MAGMA_ZGESV_GPU( magma_int_t *n, magma_int_t *nrhs, devptr_t *dA, magma_int_t *lda, magma_int_t *ipiv, devptr_t *dB, magma_int_t *ldb, magma_int_t *info)
{
  cuDoubleComplex *d_a = (cuDoubleComplex *)(*dA);
  cuDoubleComplex *d_b = (cuDoubleComplex *)(*dB);
  magma_zgesv_gpu( *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info);
}

void MAGMA_ZLARFB_GPU( char *side, char *trans, char *direct, char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *dv, magma_int_t *ldv, devptr_t *dt, magma_int_t *ldt, devptr_t *dc, magma_int_t *ldc, devptr_t *dwork, magma_int_t *ldwork)
{ 
    cuDoubleComplex *d_v = (cuDoubleComplex *)(*dv);
    cuDoubleComplex *d_t = (cuDoubleComplex *)(*dt);
    cuDoubleComplex *d_c = (cuDoubleComplex *)(*dc);
    cuDoubleComplex *d_w = (cuDoubleComplex *)(*dwork);
    magma_zlarfb_gpu( side[0], trans[0], direct[0], storev[0], *m, *n, *k, d_v, *ldv, d_t, *ldt, d_c, *ldc, d_w, *ldwork); 
}

void MAGMA_ZPOTRF_GPU( char *uplo,  magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *info)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    magma_zpotrf_gpu( uplo[0], *n, d_a, *lda, info); 
}

void MAGMA_ZPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    cuDoubleComplex *d_b = (cuDoubleComplex *)(*b);
    magma_zpotrs_gpu( uplo[0], *n, *nrhs, d_a, *lda, d_b, *ldb, info); 
}

#ifdef __cplusplus
}
#endif
