/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d

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

#undef COMPLEX
#define REAL

#ifdef ADD_

    #define MAGMA_DGEBRD magma_dgebrd_
    #define MAGMA_DGEHRD2 magma_dgehrd2_
    #define MAGMA_DGELQF magma_dgelqf_
    #define MAGMA_DGEQLF magma_dgeqlf_
    #define MAGMA_DGEQRF magma_dgeqrf_
    #define MAGMA_DGETRF magma_dgetrf_
    #define MAGMA_DLAHR2 magma_dlahr2_
    #define MAGMA_DLAHRU magma_dlahru_
    #define MAGMA_DPOTRF magma_dpotrf_
    #define MAGMA_DSYTRD magma_dsytrd_
    
    #define MAGMA_DORMQR_GPU magma_dormqr_gpu_
    #define MAGMA_DGEQRF_GPU  magma_dgeqrf_gpu_
    #define MAGMA_DGEQRF2_GPU magma_dgeqrf2_gpu_
    #define MAGMA_DGEQRS_GPU magma_dgeqrs_gpu_
    #define MAGMA_DGETRF_GPU magma_dgetrf_gpu_
    #define MAGMA_DGETRS_GPU magma_dgetrs_gpu_
    #define MAGMA_DGESV_GPU  magma_dgesv_gpu_
    #define MAGMA_DLABRD_GPU magma_dlabrd_gpu_
    #define MAGMA_DLARFB_GPU magma_dlarfb_gpu_
    #define MAGMA_DPOTRF_GPU magma_dpotrf_gpu_
    #define MAGMA_DPOTRS_GPU magma_dpotrs_gpu_

#elif defined (NOCHANGE)

    #define MAGMA_DGEBRD magma_dgebrd
    #define MAGMA_DGEHRD2 magma_dgehrd2
    #define MAGMA_DGELQF magma_dgelqf
    #define MAGMA_DGEQLF magma_dgeqlf
    #define MAGMA_DGEQRF magma_dgeqrf
    #define MAGMA_DGETRF magma_dgetrf
    #define MAGMA_DLAHR2 magma_dlahr2
    #define MAGMA_DLAHRU magma_dlahru
    #define MAGMA_DPOTRF magma_dpotrf
    #define MAGMA_DSYTRD magma_dsytrd
    
    #define MAGMA_DORMQR_GPU magma_dormqr_gpu
    #define MAGMA_DGEQRF_GPU  magma_dgeqrf_gpu
    #define MAGMA_DGEQRF2_GPU magma_dgeqrf2_gpu
    #define MAGMA_DGEQRS_GPU magma_dgeqrs_gpu
    #define MAGMA_DGETRF_GPU magma_dgetrf_gpu
    #define MAGMA_DGETRS_GPU magma_dgetrs_gpu
    #define MAGMA_DGESV_GPU  magma_dgesv_gpu
    #define MAGMA_DLABRD_GPU magma_dlabrd_gpu
    #define MAGMA_DLARFB_GPU magma_dlarfb_gpu
    #define MAGMA_DPOTRF_GPU magma_dpotrf_gpu
    #define MAGMA_DPOTRS_GPU magma_dpotrs_gpu

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
void MAGMA_DGEBRD( magma_int_t *m, magma_int_t *n, double *A, magma_int_t *lda, double *d, double *e, double *tauq,  double *taup, double *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_dgebrd( *m, *n, A, *lda, d, e, tauq,  taup, work, *lwork, info); 
}

void MAGMA_DGEHRD2( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, double *A, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_dgehrd2( *n, *ilo, *ihi, A, *lda, tau, work, lwork, info); 
}

void MAGMA_DGELQF( magma_int_t *m, magma_int_t *n, double *A, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_dgelqf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_DGEQLF( magma_int_t *m, magma_int_t *n, double *A, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_dgeqlf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_DGEQRF( magma_int_t *m, magma_int_t *n, double *A, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_dgeqrf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_DGETRF( magma_int_t *m, magma_int_t *n, double *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    magma_dgetrf( *m, *n, A, *lda, ipiv, info); 
}

void MAGMA_DLABRD_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                       double *a, magma_int_t *lda, devptr_t *da, magma_int_t *ldda, 
                       double *d, double *e, double *tauq, double *taup, 
                       double *x, magma_int_t *ldx, devptr_t *dx, magma_int_t *lddx,
                       double *y, magma_int_t *ldy, devptr_t *dy, magma_int_t *lddy)
{ 
    double *d_a = (double *)(*da);
    double *d_x = (double *)(*dx);
    double *d_y = (double *)(*dy);
    magma_dlabrd_gpu( *m, *n, *nb, 
                      a, *lda, d_a, *ldda, 
                      d, e, tauq, taup, 
                      x, *ldx, d_x, *lddx,
                      y, *ldy, d_y, *lddy); 
}

void MAGMA_DLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, devptr_t *da, devptr_t *dv, double *a, magma_int_t *lda, double *tau, double *t, magma_int_t *ldt, double *y, magma_int_t *ldy)
{ 
    double *d_a = (double *)(*da);
    double *d_v = (double *)(*dv);
    magma_dlahr2( *m, *n, *nb, d_a, d_v, a, *lda, tau, t, *ldt, y, *ldy); 
}

void MAGMA_DLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, double *a, magma_int_t *lda, devptr_t *da, double *y, double *v, double *t, devptr_t *dwork)
{ 
    double *d_a = (double *)(*da);
    double *d_w = (double *)(*dwork);
    magma_dlahru( *m, *n, *nb, a, *lda, d_a, y, v, t, d_w); 
}

void MAGMA_DPOTRF( char *uplo, magma_int_t *n, double *A, magma_int_t *lda, magma_int_t *info)
{ 
    magma_dpotrf( uplo[0], *n, A, *lda, info); 
}

void MAGMA_DSYTRD( char *uplo, magma_int_t *n, double *A, magma_int_t *lda, double *d, double *e, double *tau, double *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_dsytrd( uplo[0], *n, A, *lda, d, e, tau, work, *lwork, info); 
}


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMA_DORMQR_GPU(char *side, char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *a, magma_int_t *lda, double *tau, devptr_t *c, magma_int_t *ldc, double *work, magma_int_t *lwork, devptr_t *td, magma_int_t *nb, magma_int_t *info)
{ 
    double *d_a  = (double *)(*a);
    double *d_c  = (double *)(*c);
    double *d_td = (double *)(*td);
    magma_dormqr_gpu(side[0], trans[0], *m, *n, *k, d_a, *lda, tau, d_c, *ldc, work, *lwork, d_td, *nb, info); 
}

void MAGMA_DGEQRF2_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, double *tau, magma_int_t *info)
{ 
    double *d_a = (double *)(*A);
    magma_dgeqrf2_gpu( *m, *n, d_a, *lda, tau, info); 
}

void MAGMA_DGEQRF_GPU(magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, double *tau, devptr_t *dwork, magma_int_t *info)
{ 
    double *d_a = (double *)(*A);
    double *d_w = (double *)(*dwork);
    magma_dgeqrf_gpu(*m, *n, d_a, *lda, tau, d_w, info); 
}

void MAGMA_DGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
		       devptr_t *A, magma_int_t *lda, double *tau, devptr_t *td, 
		       devptr_t *c, magma_int_t *ldc, double *work, magma_int_t *lwork, magma_int_t *info)
{ 
    double *d_a  = (double *)(*A);
    double *d_c  = (double *)(*c);
    double *d_td = (double *)(*td);
    magma_dgeqrs_gpu( *m, *n, *nrhs, d_a, *lda, tau, d_td, d_c, *ldc, work, *lwork, info); 
}

void MAGMA_DGETRF_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    double *d_a = (double *)(*A);
    magma_dgetrf_gpu( *m, *n, d_a, *lda, ipiv, info); 
}

void MAGMA_DGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
  double *d_a = (double *)(uintptr_t)(*A);
  double *d_b = (double *)(uintptr_t)(*b);
  magma_dgetrs_gpu( trans[0], *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info); 
}

void MAGMA_DGESV_GPU( magma_int_t *n, magma_int_t *nrhs, devptr_t *dA, magma_int_t *lda, magma_int_t *ipiv, devptr_t *dB, magma_int_t *ldb, magma_int_t *info)
{
  double *d_a = (double *)(*dA);
  double *d_b = (double *)(*dB);
  magma_dgesv_gpu( *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info);
}

void MAGMA_DLARFB_GPU( char *side, char *trans, char *direct, char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *dv, magma_int_t *ldv, devptr_t *dt, magma_int_t *ldt, devptr_t *dc, magma_int_t *ldc, devptr_t *dwork, magma_int_t *ldwork)
{ 
    double *d_v = (double *)(*dv);
    double *d_t = (double *)(*dt);
    double *d_c = (double *)(*dc);
    double *d_w = (double *)(*dwork);
    magma_dlarfb_gpu( side[0], trans[0], direct[0], storev[0], *m, *n, *k, d_v, *ldv, d_t, *ldt, d_c, *ldc, d_w, *ldwork); 
}

void MAGMA_DPOTRF_GPU( char *uplo,  magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *info)
{ 
    double *d_a = (double *)(*A);
    magma_dpotrf_gpu( uplo[0], *n, d_a, *lda, info); 
}

void MAGMA_DPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
    double *d_a = (double *)(*A);
    double *d_b = (double *)(*b);
    magma_dpotrs_gpu( uplo[0], *n, *nrhs, d_a, *lda, d_b, *ldb, info); 
}

#ifdef __cplusplus
}
#endif
