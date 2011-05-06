/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated s

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

    #define MAGMA_SGEBRD magma_sgebrd_
    #define MAGMA_SGEHRD2 magma_sgehrd2_
    #define MAGMA_SGELQF magma_sgelqf_
    #define MAGMA_SGEQLF magma_sgeqlf_
    #define MAGMA_SGEQRF magma_sgeqrf_
    #define MAGMA_SGETRF magma_sgetrf_
    #define MAGMA_SLAHR2 magma_slahr2_
    #define MAGMA_SLAHRU magma_slahru_
    #define MAGMA_SPOTRF magma_spotrf_
    #define MAGMA_SSYTRD magma_ssytrd_
    
    #define MAGMA_SORMQR_GPU magma_sormqr_gpu_
    #define MAGMA_SGEQRF_GPU  magma_sgeqrf_gpu_
    #define MAGMA_SGEQRF2_GPU magma_sgeqrf2_gpu_
    #define MAGMA_SGEQRS_GPU magma_sgeqrs_gpu_
    #define MAGMA_SGETRF_GPU magma_sgetrf_gpu_
    #define MAGMA_SGETRS_GPU magma_sgetrs_gpu_
    #define MAGMA_SGESV_GPU  magma_sgesv_gpu_
    #define MAGMA_SLABRD_GPU magma_slabrd_gpu_
    #define MAGMA_SLARFB_GPU magma_slarfb_gpu_
    #define MAGMA_SPOTRF_GPU magma_spotrf_gpu_
    #define MAGMA_SPOTRS_GPU magma_spotrs_gpu_

#elif defined (NOCHANGE)

    #define MAGMA_SGEBRD magma_sgebrd
    #define MAGMA_SGEHRD2 magma_sgehrd2
    #define MAGMA_SGELQF magma_sgelqf
    #define MAGMA_SGEQLF magma_sgeqlf
    #define MAGMA_SGEQRF magma_sgeqrf
    #define MAGMA_SGETRF magma_sgetrf
    #define MAGMA_SLAHR2 magma_slahr2
    #define MAGMA_SLAHRU magma_slahru
    #define MAGMA_SPOTRF magma_spotrf
    #define MAGMA_SSYTRD magma_ssytrd
    
    #define MAGMA_SORMQR_GPU magma_sormqr_gpu
    #define MAGMA_SGEQRF_GPU  magma_sgeqrf_gpu
    #define MAGMA_SGEQRF2_GPU magma_sgeqrf2_gpu
    #define MAGMA_SGEQRS_GPU magma_sgeqrs_gpu
    #define MAGMA_SGETRF_GPU magma_sgetrf_gpu
    #define MAGMA_SGETRS_GPU magma_sgetrs_gpu
    #define MAGMA_SGESV_GPU  magma_sgesv_gpu
    #define MAGMA_SLABRD_GPU magma_slabrd_gpu
    #define MAGMA_SLARFB_GPU magma_slarfb_gpu
    #define MAGMA_SPOTRF_GPU magma_spotrf_gpu
    #define MAGMA_SPOTRS_GPU magma_spotrs_gpu

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
void MAGMA_SGEBRD( magma_int_t *m, magma_int_t *n, float *A, magma_int_t *lda, float *d, float *e, float *tauq,  float *taup, float *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_sgebrd( *m, *n, A, *lda, d, e, tauq,  taup, work, *lwork, info); 
}

void MAGMA_SGEHRD2( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, float *A, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_sgehrd2( *n, *ilo, *ihi, A, *lda, tau, work, lwork, info); 
}

void MAGMA_SGELQF( magma_int_t *m, magma_int_t *n, float *A, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_sgelqf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_SGEQLF( magma_int_t *m, magma_int_t *n, float *A, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_sgeqlf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_SGEQRF( magma_int_t *m, magma_int_t *n, float *A, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_sgeqrf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_SGETRF( magma_int_t *m, magma_int_t *n, float *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    magma_sgetrf( *m, *n, A, *lda, ipiv, info); 
}

void MAGMA_SLABRD_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                       float *a, magma_int_t *lda, devptr_t *da, magma_int_t *ldda, 
                       float *d, float *e, float *tauq, float *taup, 
                       float *x, magma_int_t *ldx, devptr_t *dx, magma_int_t *lddx,
                       float *y, magma_int_t *ldy, devptr_t *dy, magma_int_t *lddy)
{ 
    float *d_a = (float *)(*da);
    float *d_x = (float *)(*dx);
    float *d_y = (float *)(*dy);
    magma_slabrd_gpu( *m, *n, *nb, 
                      a, *lda, d_a, *ldda, 
                      d, e, tauq, taup, 
                      x, *ldx, d_x, *lddx,
                      y, *ldy, d_y, *lddy); 
}

void MAGMA_SLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, devptr_t *da, devptr_t *dv, float *a, magma_int_t *lda, float *tau, float *t, magma_int_t *ldt, float *y, magma_int_t *ldy)
{ 
    float *d_a = (float *)(*da);
    float *d_v = (float *)(*dv);
    magma_slahr2( *m, *n, *nb, d_a, d_v, a, *lda, tau, t, *ldt, y, *ldy); 
}

void MAGMA_SLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, float *a, magma_int_t *lda, devptr_t *da, float *y, float *v, float *t, devptr_t *dwork)
{ 
    float *d_a = (float *)(*da);
    float *d_w = (float *)(*dwork);
    magma_slahru( *m, *n, *nb, a, *lda, d_a, y, v, t, d_w); 
}

void MAGMA_SPOTRF( char *uplo, magma_int_t *n, float *A, magma_int_t *lda, magma_int_t *info)
{ 
    magma_spotrf( uplo[0], *n, A, *lda, info); 
}

void MAGMA_SSYTRD( char *uplo, magma_int_t *n, float *A, magma_int_t *lda, float *d, float *e, float *tau, float *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_ssytrd( uplo[0], *n, A, *lda, d, e, tau, work, *lwork, info); 
}


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMA_SORMQR_GPU(char *side, char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *a, magma_int_t *lda, float *tau, devptr_t *c, magma_int_t *ldc, float *work, magma_int_t *lwork, devptr_t *td, magma_int_t *nb, magma_int_t *info)
{ 
    float *d_a  = (float *)(*a);
    float *d_c  = (float *)(*c);
    float *d_td = (float *)(*td);
    magma_sormqr_gpu(side[0], trans[0], *m, *n, *k, d_a, *lda, tau, d_c, *ldc, work, *lwork, d_td, *nb, info); 
}

void MAGMA_SGEQRF2_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, float *tau, magma_int_t *info)
{ 
    float *d_a = (float *)(*A);
    magma_sgeqrf2_gpu( *m, *n, d_a, *lda, tau, info); 
}

void MAGMA_SGEQRF_GPU(magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, float *tau, devptr_t *dwork, magma_int_t *info)
{ 
    float *d_a = (float *)(*A);
    float *d_w = (float *)(*dwork);
    magma_sgeqrf_gpu(*m, *n, d_a, *lda, tau, d_w, info); 
}

void MAGMA_SGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
		       devptr_t *A, magma_int_t *lda, float *tau, devptr_t *td, 
		       devptr_t *c, magma_int_t *ldc, float *work, magma_int_t *lwork, magma_int_t *info)
{ 
    float *d_a  = (float *)(*A);
    float *d_c  = (float *)(*c);
    float *d_td = (float *)(*td);
    magma_sgeqrs_gpu( *m, *n, *nrhs, d_a, *lda, tau, d_td, d_c, *ldc, work, *lwork, info); 
}

void MAGMA_SGETRF_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    float *d_a = (float *)(*A);
    magma_sgetrf_gpu( *m, *n, d_a, *lda, ipiv, info); 
}

void MAGMA_SGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
  float *d_a = (float *)(uintptr_t)(*A);
  float *d_b = (float *)(uintptr_t)(*b);
  magma_sgetrs_gpu( trans[0], *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info); 
}

void MAGMA_SGESV_GPU( magma_int_t *n, magma_int_t *nrhs, devptr_t *dA, magma_int_t *lda, magma_int_t *ipiv, devptr_t *dB, magma_int_t *ldb, magma_int_t *info)
{
  float *d_a = (float *)(*dA);
  float *d_b = (float *)(*dB);
  magma_sgesv_gpu( *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info);
}

void MAGMA_SLARFB_GPU( char *side, char *trans, char *direct, char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *dv, magma_int_t *ldv, devptr_t *dt, magma_int_t *ldt, devptr_t *dc, magma_int_t *ldc, devptr_t *dwork, magma_int_t *ldwork)
{ 
    float *d_v = (float *)(*dv);
    float *d_t = (float *)(*dt);
    float *d_c = (float *)(*dc);
    float *d_w = (float *)(*dwork);
    magma_slarfb_gpu( side[0], trans[0], direct[0], storev[0], *m, *n, *k, d_v, *ldv, d_t, *ldt, d_c, *ldc, d_w, *ldwork); 
}

void MAGMA_SPOTRF_GPU( char *uplo,  magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *info)
{ 
    float *d_a = (float *)(*A);
    magma_spotrf_gpu( uplo[0], *n, d_a, *lda, info); 
}

void MAGMA_SPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
    float *d_a = (float *)(*A);
    float *d_b = (float *)(*b);
    magma_spotrs_gpu( uplo[0], *n, *nrhs, d_a, *lda, d_b, *ldb, info); 
}

#ifdef __cplusplus
}
#endif
