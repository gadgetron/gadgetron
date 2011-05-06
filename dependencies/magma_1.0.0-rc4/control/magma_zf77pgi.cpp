/**
 *
 * @file magma_zf77.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.1.0
 * @author Mathieu Faverge
 * @date 2009-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <stdlib.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include "magma.h"

#undef REAL
#define COMPLEX

#ifdef ADD_

    #define MAGMA_ZGEBRD magma_zgebrd_
    #define MAGMA_ZGEHRD magma_zgehrd_
    #define MAGMA_ZGELQF magma_zgelqf_
    #define MAGMA_ZGEQLF magma_zgeqlf_
    #define MAGMA_ZGEQRF magma_zgeqrf_
    #define MAGMA_ZGETRF magma_zgetrf_
    #define MAGMA_ZLABRD magma_zlabrd_
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
    #define MAGMA_ZLARFB_GPU magma_zlarfb_gpu_
    #define MAGMA_ZPOTRF_GPU magma_zpotrf_gpu_
    #define MAGMA_ZPOTRS_GPU magma_zpotrs_gpu_

#elif defined (NOCHANGE)

    #define MAGMA_ZGEBRD magma_zgebrd
    #define MAGMA_ZGEHRD magma_zgehrd
    #define MAGMA_ZGELQF magma_zgelqf
    #define MAGMA_ZGEQLF magma_zgeqlf
    #define MAGMA_ZGEQRF magma_zgeqrf
    #define MAGMA_ZGETRF magma_zgetrf
    #define MAGMA_ZLABRD magma_zlabrd
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
void MAGMA_ZGEBRD( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double *d, double *e, double2 *tauq,  double2 *taup, double2 *work, magma_int_t *lwork, double2 *da, magma_int_t *info)
{ magma_zgebrd( *m, *n, A, *lda, d, e, tauq,  taup, work, *lwork, da, info); }

void MAGMA_ZGEHRD( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, double2 *da, magma_int_t *info)
{ magma_zgehrd( *n, *ilo, *ihi, A, *lda, tau, work, lwork, da, info); }

void MAGMA_ZGELQF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ magma_zgelqf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_ZGEQLF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ magma_zgeqlf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_ZGEQRF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ magma_zgeqrf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_ZGETRF( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ magma_zgetrf( *m, *n, A, *lda, ipiv, info); }

void MAGMA_ZLABRD( magma_int_t *m, magma_int_t *n, magma_int_t *nb, double2 *a, magma_int_t *lda, double *d, double *e, double2 *tauq, double2 *taup, double2 *x, magma_int_t *ldx, double2 *y, magma_int_t *ldy, double2 *da, magma_int_t *ldda, double2 *dx, magma_int_t *lddx, double2 *dy, magma_int_t *lddy)
{ magma_zlabrd( *m, *n, *nb, a, *lda, d, e, tauq, taup, x, *ldx, y, *ldy, da, *ldda, dx, *lddx, dy, *lddy); }

void MAGMA_ZLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, double2 *da, double2 *dv, double2 *a, magma_int_t *lda, double2 *tau, double2 *t, magma_int_t *ldt, double2 *y, magma_int_t *ldy)
{ magma_zlahr2( *m, *n, *nb, da, dv, a, *lda, tau, t, *ldt, y, *ldy); }

void MAGMA_ZLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, double2 *a, magma_int_t *lda, double2 *da, double2 *y, double2 *v, double2 *t, double2 *dwork)
{ magma_zlahru( *m, *n, *nb, a, *lda, da, y, v, t, dwork); }

void MAGMA_ZPOTRF( char *uplo, magma_int_t *n, double2 *A, magma_int_t *lda, magma_int_t *info)
{ magma_zpotrf( uplo[0], *n, A, *lda, info); }

void MAGMA_ZHETRD( char *uplo, magma_int_t *n, double2 *A, magma_int_t *lda, double *d, double *e, double2 *tau, double2 *work, magma_int_t *lwork, double2 *da, magma_int_t *info)
{ magma_zhetrd( uplo[0], *n, A, *lda, d, e, tau, work, *lwork, da, info); }


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMA_ZUNMQR_GPU(char *side, char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, double2 *a, magma_int_t *lda, double2 *tau, double2 *c, magma_int_t *ldc, double2 *work, magma_int_t *lwork, double2 *td, magma_int_t *nb, magma_int_t *info)
{ magma_zunmqr_gpu(side[0], trans[0], *m, *n, *k, a, *lda, tau, c, *ldc, work, *lwork, td, *nb, info); }

void MAGMA_ZGEQRF2_GPU( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, magma_int_t *info)
{ magma_zgeqrf2_gpu( *m, *n, A, *lda, tau, info); }

void MAGMA_ZGEQRF_GPU(magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, double2 *tau, double2 *dwork, magma_int_t *info)
{ magma_zgeqrf_gpu(*m, *n, A, *lda, tau, dwork, info); }

void MAGMA_ZGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, double2 *A, magma_int_t *lda, double2 *tau, double2 *td, double2 *c, magma_int_t *ldc, double2 *work, magma_int_t *lwork, magma_int_t *info)
{ magma_zgeqrs_gpu( *m, *n, *nrhs, A, *lda, tau, td, c, *ldc, work, *lwork, info); }

void MAGMA_ZGETRF_GPU( magma_int_t *m, magma_int_t *n, double2 *A, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info)
{ magma_zgetrf_gpu( *m, *n, A, *lda, ipiv, info); }

void MAGMA_ZGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, double2 *A, magma_int_t *lda, magma_int_t *ipiv, double2 *b, magma_int_t *ldb, magma_int_t *info)
{ magma_zgetrs_gpu( trans[0], *n, *nrhs, A, *lda, ipiv, b, *ldb, info ); }

void MAGMA_ZLARFB_GPU( char *side, char *trans, char *direct, char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, double2 *dv, magma_int_t *ldv, double2 *dt, magma_int_t *ldt, double2 *dc, magma_int_t *ldc, double2 *dowrk, magma_int_t *ldwork)
{ magma_zlarfb_gpu( side[0], trans[0], direct[0], storev[0], *m, *n, *k, dv, *ldv, dt, *ldt, dc, *ldc, dowrk, *ldwork); }

void MAGMA_ZPOTRF_GPU( char *uplo,  magma_int_t *n, double2 *A, magma_int_t *lda, magma_int_t *info)
{ magma_zpotrf_gpu( uplo[0], *n, A, *lda, info); }

void MAGMA_ZPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, double2 *A, magma_int_t *lda, double2 *b, magma_int_t *ldb, magma_int_t *info)
{ magma_zpotrs_gpu( uplo[0], *n, *nrhs, A, *lda, b, *ldb, info); }

#ifdef __cplusplus
}
#endif
