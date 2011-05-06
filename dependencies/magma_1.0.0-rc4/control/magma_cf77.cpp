/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated c

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

    #define MAGMA_CGEBRD magma_cgebrd_
    #define MAGMA_CGEHRD2 magma_cgehrd2_
    #define MAGMA_CGELQF magma_cgelqf_
    #define MAGMA_CGEQLF magma_cgeqlf_
    #define MAGMA_CGEQRF magma_cgeqrf_
    #define MAGMA_CGETRF magma_cgetrf_
    #define MAGMA_CLAHR2 magma_clahr2_
    #define MAGMA_CLAHRU magma_clahru_
    #define MAGMA_CPOTRF magma_cpotrf_
    #define MAGMA_CHETRD magma_chetrd_
    
    #define MAGMA_CUNMQR_GPU magma_cunmqr_gpu_
    #define MAGMA_CGEQRF_GPU  magma_cgeqrf_gpu_
    #define MAGMA_CGEQRF2_GPU magma_cgeqrf2_gpu_
    #define MAGMA_CGEQRS_GPU magma_cgeqrs_gpu_
    #define MAGMA_CGETRF_GPU magma_cgetrf_gpu_
    #define MAGMA_CGETRS_GPU magma_cgetrs_gpu_
    #define MAGMA_CGESV_GPU  magma_cgesv_gpu_
    #define MAGMA_CLABRD_GPU magma_clabrd_gpu_
    #define MAGMA_CLARFB_GPU magma_clarfb_gpu_
    #define MAGMA_CPOTRF_GPU magma_cpotrf_gpu_
    #define MAGMA_CPOTRS_GPU magma_cpotrs_gpu_

#elif defined (NOCHANGE)

    #define MAGMA_CGEBRD magma_cgebrd
    #define MAGMA_CGEHRD2 magma_cgehrd2
    #define MAGMA_CGELQF magma_cgelqf
    #define MAGMA_CGEQLF magma_cgeqlf
    #define MAGMA_CGEQRF magma_cgeqrf
    #define MAGMA_CGETRF magma_cgetrf
    #define MAGMA_CLAHR2 magma_clahr2
    #define MAGMA_CLAHRU magma_clahru
    #define MAGMA_CPOTRF magma_cpotrf
    #define MAGMA_CHETRD magma_chetrd
    
    #define MAGMA_CUNMQR_GPU magma_cunmqr_gpu
    #define MAGMA_CGEQRF_GPU  magma_cgeqrf_gpu
    #define MAGMA_CGEQRF2_GPU magma_cgeqrf2_gpu
    #define MAGMA_CGEQRS_GPU magma_cgeqrs_gpu
    #define MAGMA_CGETRF_GPU magma_cgetrf_gpu
    #define MAGMA_CGETRS_GPU magma_cgetrs_gpu
    #define MAGMA_CGESV_GPU  magma_cgesv_gpu
    #define MAGMA_CLABRD_GPU magma_clabrd_gpu
    #define MAGMA_CLARFB_GPU magma_clarfb_gpu
    #define MAGMA_CPOTRF_GPU magma_cpotrf_gpu
    #define MAGMA_CPOTRS_GPU magma_cpotrs_gpu

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
void MAGMA_CGEBRD( magma_int_t *m, magma_int_t *n, float2 *A, magma_int_t *lda, float *d, float *e, float2 *tauq,  float2 *taup, float2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_cgebrd( *m, *n, A, *lda, d, e, tauq,  taup, work, *lwork, info); 
}

void MAGMA_CGEHRD2( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, float2 *A, magma_int_t *lda, float2 *tau, float2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_cgehrd2( *n, *ilo, *ihi, A, *lda, tau, work, lwork, info); 
}

void MAGMA_CGELQF( magma_int_t *m, magma_int_t *n, float2 *A, magma_int_t *lda, float2 *tau, float2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_cgelqf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_CGEQLF( magma_int_t *m, magma_int_t *n, float2 *A, magma_int_t *lda, float2 *tau, float2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_cgeqlf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_CGEQRF( magma_int_t *m, magma_int_t *n, float2 *A, magma_int_t *lda, float2 *tau, float2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_cgeqrf( *m, *n, A, *lda, tau, work, *lwork, info); 
}

void MAGMA_CGETRF( magma_int_t *m, magma_int_t *n, float2 *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    magma_cgetrf( *m, *n, A, *lda, ipiv, info); 
}

void MAGMA_CLABRD_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                       float2 *a, magma_int_t *lda, devptr_t *da, magma_int_t *ldda, 
                       float *d, float *e, float2 *tauq, float2 *taup, 
                       float2 *x, magma_int_t *ldx, devptr_t *dx, magma_int_t *lddx,
                       float2 *y, magma_int_t *ldy, devptr_t *dy, magma_int_t *lddy)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*da);
    cuFloatComplex *d_x = (cuFloatComplex *)(*dx);
    cuFloatComplex *d_y = (cuFloatComplex *)(*dy);
    magma_clabrd_gpu( *m, *n, *nb, 
                      a, *lda, d_a, *ldda, 
                      d, e, tauq, taup, 
                      x, *ldx, d_x, *lddx,
                      y, *ldy, d_y, *lddy); 
}

void MAGMA_CLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, devptr_t *da, devptr_t *dv, float2 *a, magma_int_t *lda, float2 *tau, float2 *t, magma_int_t *ldt, float2 *y, magma_int_t *ldy)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*da);
    cuFloatComplex *d_v = (cuFloatComplex *)(*dv);
    magma_clahr2( *m, *n, *nb, d_a, d_v, a, *lda, tau, t, *ldt, y, *ldy); 
}

void MAGMA_CLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, float2 *a, magma_int_t *lda, devptr_t *da, float2 *y, float2 *v, float2 *t, devptr_t *dwork)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*da);
    cuFloatComplex *d_w = (cuFloatComplex *)(*dwork);
    magma_clahru( *m, *n, *nb, a, *lda, d_a, y, v, t, d_w); 
}

void MAGMA_CPOTRF( char *uplo, magma_int_t *n, float2 *A, magma_int_t *lda, magma_int_t *info)
{ 
    magma_cpotrf( uplo[0], *n, A, *lda, info); 
}

void MAGMA_CHETRD( char *uplo, magma_int_t *n, float2 *A, magma_int_t *lda, float *d, float *e, float2 *tau, float2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    magma_chetrd( uplo[0], *n, A, *lda, d, e, tau, work, *lwork, info); 
}


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMA_CUNMQR_GPU(char *side, char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *a, magma_int_t *lda, float2 *tau, devptr_t *c, magma_int_t *ldc, float2 *work, magma_int_t *lwork, devptr_t *td, magma_int_t *nb, magma_int_t *info)
{ 
    cuFloatComplex *d_a  = (cuFloatComplex *)(*a);
    cuFloatComplex *d_c  = (cuFloatComplex *)(*c);
    cuFloatComplex *d_td = (cuFloatComplex *)(*td);
    magma_cunmqr_gpu(side[0], trans[0], *m, *n, *k, d_a, *lda, tau, d_c, *ldc, work, *lwork, d_td, *nb, info); 
}

void MAGMA_CGEQRF2_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, float2 *tau, magma_int_t *info)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*A);
    magma_cgeqrf2_gpu( *m, *n, d_a, *lda, tau, info); 
}

void MAGMA_CGEQRF_GPU(magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, float2 *tau, devptr_t *dwork, magma_int_t *info)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*A);
    cuFloatComplex *d_w = (cuFloatComplex *)(*dwork);
    magma_cgeqrf_gpu(*m, *n, d_a, *lda, tau, d_w, info); 
}

void MAGMA_CGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
		       devptr_t *A, magma_int_t *lda, float2 *tau, devptr_t *td, 
		       devptr_t *c, magma_int_t *ldc, float2 *work, magma_int_t *lwork, magma_int_t *info)
{ 
    cuFloatComplex *d_a  = (cuFloatComplex *)(*A);
    cuFloatComplex *d_c  = (cuFloatComplex *)(*c);
    cuFloatComplex *d_td = (cuFloatComplex *)(*td);
    magma_cgeqrs_gpu( *m, *n, *nrhs, d_a, *lda, tau, d_td, d_c, *ldc, work, *lwork, info); 
}

void MAGMA_CGETRF_GPU( magma_int_t *m, magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*A);
    magma_cgetrf_gpu( *m, *n, d_a, *lda, ipiv, info); 
}

void MAGMA_CGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, magma_int_t *ipiv, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
  cuFloatComplex *d_a = (cuFloatComplex *)(uintptr_t)(*A);
  cuFloatComplex *d_b = (cuFloatComplex *)(uintptr_t)(*b);
  magma_cgetrs_gpu( trans[0], *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info); 
}

void MAGMA_CGESV_GPU( magma_int_t *n, magma_int_t *nrhs, devptr_t *dA, magma_int_t *lda, magma_int_t *ipiv, devptr_t *dB, magma_int_t *ldb, magma_int_t *info)
{
  cuFloatComplex *d_a = (cuFloatComplex *)(*dA);
  cuFloatComplex *d_b = (cuFloatComplex *)(*dB);
  magma_cgesv_gpu( *n, *nrhs, d_a, *lda, ipiv, d_b, *ldb, info);
}

void MAGMA_CLARFB_GPU( char *side, char *trans, char *direct, char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, devptr_t *dv, magma_int_t *ldv, devptr_t *dt, magma_int_t *ldt, devptr_t *dc, magma_int_t *ldc, devptr_t *dwork, magma_int_t *ldwork)
{ 
    cuFloatComplex *d_v = (cuFloatComplex *)(*dv);
    cuFloatComplex *d_t = (cuFloatComplex *)(*dt);
    cuFloatComplex *d_c = (cuFloatComplex *)(*dc);
    cuFloatComplex *d_w = (cuFloatComplex *)(*dwork);
    magma_clarfb_gpu( side[0], trans[0], direct[0], storev[0], *m, *n, *k, d_v, *ldv, d_t, *ldt, d_c, *ldc, d_w, *ldwork); 
}

void MAGMA_CPOTRF_GPU( char *uplo,  magma_int_t *n, devptr_t *A, magma_int_t *lda, magma_int_t *info)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*A);
    magma_cpotrf_gpu( uplo[0], *n, d_a, *lda, info); 
}

void MAGMA_CPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, devptr_t *A, magma_int_t *lda, devptr_t *b, magma_int_t *ldb, magma_int_t *info)
{ 
    cuFloatComplex *d_a = (cuFloatComplex *)(*A);
    cuFloatComplex *d_b = (cuFloatComplex *)(*b);
    magma_cpotrs_gpu( uplo[0], *n, *nrhs, d_a, *lda, d_b, *ldb, info); 
}

#ifdef __cplusplus
}
#endif
