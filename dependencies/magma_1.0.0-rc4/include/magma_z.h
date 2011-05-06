/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMA_Z_H_
#define _MAGMA_Z_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_zgebrd( magma_int_t m, magma_int_t n, cuDoubleComplex *A, 
			  magma_int_t lda, double *d, double *e,
			  cuDoubleComplex *tauq,  cuDoubleComplex *taup, 
			  cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  cuDoubleComplex *A, magma_int_t lda, cuDoubleComplex *tau, 
			  cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_zgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  cuDoubleComplex *A, magma_int_t lda, cuDoubleComplex *tau,
			  cuDoubleComplex *work, magma_int_t lwork,
			  cuDoubleComplex *d_T, magma_int_t *info);
magma_int_t magma_zgelqf( magma_int_t m, magma_int_t n, 
                          cuDoubleComplex *A,    magma_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgeqlf( magma_int_t m, magma_int_t n, 
                          cuDoubleComplex *A,    magma_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgeqrf( magma_int_t m, magma_int_t n, cuDoubleComplex *A, 
			  magma_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
			  magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgetrf( magma_int_t m, magma_int_t n, cuDoubleComplex *A, 
			  magma_int_t lda, magma_int_t *ipiv, 
			  magma_int_t *info);
magma_int_t magma_zlatrd( char uplo, magma_int_t n, magma_int_t nb, cuDoubleComplex *a, 
                          magma_int_t lda, double *e, cuDoubleComplex *tau, 
			  cuDoubleComplex *w, magma_int_t ldw,
                          cuDoubleComplex *da, magma_int_t ldda, 
			  cuDoubleComplex *dw, magma_int_t lddw);
magma_int_t magma_zlahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  cuDoubleComplex *da, cuDoubleComplex *dv, cuDoubleComplex *a, 
			  magma_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *t, 
			  magma_int_t ldt, cuDoubleComplex *y, magma_int_t ldy);
magma_int_t magma_zlahru( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  cuDoubleComplex *a, magma_int_t lda, 
			  cuDoubleComplex *da, cuDoubleComplex *y, 
			  cuDoubleComplex *v, cuDoubleComplex *t, 
			  cuDoubleComplex *dwork);
magma_int_t magma_zhebbd( char uplo, magma_int_t n, cuDoubleComplex *a, magma_int_t lda,
			  cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_zpotrf( char uplo, magma_int_t n, cuDoubleComplex *A, 
			  magma_int_t lda, magma_int_t *info);
magma_int_t magma_zhetrd( char uplo, magma_int_t n, cuDoubleComplex *A, 
			  magma_int_t lda, double *d, double *e, 
			  cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t lwork, 
			  magma_int_t *info);
magma_int_t magma_zungqr( magma_int_t m, magma_int_t n, magma_int_t k,
			  cuDoubleComplex *a, magma_int_t lda,
			  cuDoubleComplex *tau, cuDoubleComplex *dwork,
			  magma_int_t nb, magma_int_t *info );
magma_int_t magma_zunmqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex *tau, 
                          cuDoubleComplex *c, magma_int_t ldc, 
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zunmtr( char side, char uplo, char trans,
			  int m, int n,
			  cuDoubleComplex *a,    int lda,
			  cuDoubleComplex *tau,
			  cuDoubleComplex *c,    int ldc,
			  cuDoubleComplex *work, int lwork,
			  int *info);
magma_int_t magma_zunghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  cuDoubleComplex *a, magma_int_t lda,
			  cuDoubleComplex *tau,
			  cuDoubleComplex *dT, magma_int_t nb,
			  magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_zgeev( char jobvl, char jobvr, magma_int_t n,
			  cuDoubleComplex *a, magma_int_t lda,
			  cuDoubleComplex *w,
			  cuDoubleComplex *vl, magma_int_t ldvl,
			  cuDoubleComplex *vr, magma_int_t ldvr,
			  cuDoubleComplex *work, magma_int_t lwork,
			  double *rwork, magma_int_t *info);
magma_int_t magma_zgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  cuDoubleComplex *a,    magma_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_int_t ldu, 
			  cuDoubleComplex *vt,   magma_int_t ldvt,
			  cuDoubleComplex *work, magma_int_t lwork,
			  double *rwork, magma_int_t *info );
magma_int_t magma_zheevd( char jobz, char uplo, magma_int_t n,
			  cuDoubleComplex *a, magma_int_t lda, double *w,
			  cuDoubleComplex *work, magma_int_t lwork,
		          double *rwork, magma_int_t lrwork,
		          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#else
magma_int_t  magma_zgeev( char jobvl, char jobvr, magma_int_t n,
			  cuDoubleComplex *a,    magma_int_t lda,
			  cuDoubleComplex *wr, cuDoubleComplex *wi,
			  cuDoubleComplex *vl,   magma_int_t ldvl,
			  cuDoubleComplex *vr,   magma_int_t ldvr,
			  cuDoubleComplex *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_zgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  cuDoubleComplex *a,    magma_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_int_t ldu, 
			  cuDoubleComplex *vt,   magma_int_t ldvt,
			  cuDoubleComplex *work, magma_int_t lwork,
			  magma_int_t *info );
magma_int_t magma_zheevd( char jobz, char uplo, magma_int_t n,
			  cuDoubleComplex *a, magma_int_t lda, double *w,
			  cuDoubleComplex *work, magma_int_t lwork,
			  magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_zgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
			      cuDoubleComplex *dA,    magma_int_t ldda, 
			      cuDoubleComplex *dB,    magma_int_t lddb, 
			      cuDoubleComplex *hwork, magma_int_t lwork, 
			      magma_int_t *info);
magma_int_t magma_zgeqrf_gpu( magma_int_t m, magma_int_t n, 
			      cuDoubleComplex *dA,  magma_int_t ldda, 
			      cuDoubleComplex *tau, cuDoubleComplex *dT, 
			      magma_int_t *info);
magma_int_t magma_zgeqrf2_gpu(magma_int_t m, magma_int_t n, 
			      cuDoubleComplex *dA,  magma_int_t ldda, 
			      cuDoubleComplex *tau, magma_int_t *info);
magma_int_t magma_zgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
			      cuDoubleComplex *dA,     magma_int_t ldda, 
			      cuDoubleComplex *tau,   cuDoubleComplex *dT,
			      cuDoubleComplex *dB,    magma_int_t lddb,
			      cuDoubleComplex *hwork, magma_int_t lhwork, 
			      magma_int_t *info);
magma_int_t magma_zgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              cuDoubleComplex *dL1, magma_int_t lddl1, 
                              cuDoubleComplex *dL,  magma_int_t lddl, 
                              cuDoubleComplex *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_zgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
			      cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zgetfl_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda,
                              cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuDoubleComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_zgetrf_gpu( magma_int_t m, magma_int_t n, 
			      cuDoubleComplex *dA, magma_int_t ldda, 
			      magma_int_t *ipiv, magma_int_t *info);
magma_int_t magma_zgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
			      cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zlabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex *da, magma_int_t ldda,
                              double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,  
                              cuDoubleComplex *x, magma_int_t ldx, cuDoubleComplex *dx, magma_int_t lddx, 
                              cuDoubleComplex *y, magma_int_t ldy, cuDoubleComplex *dy, magma_int_t lddy);
magma_int_t magma_zlarfb_gpu( char side, char trans, char direct, char storev, 
			      magma_int_t m, magma_int_t n, magma_int_t k,
			      cuDoubleComplex *dv, magma_int_t ldv, cuDoubleComplex *dt,    magma_int_t ldt, 
			      cuDoubleComplex *dc, magma_int_t ldc, cuDoubleComplex *dowrk, magma_int_t ldwork );
magma_int_t magma_zposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
			      cuDoubleComplex *dA, magma_int_t ldda, 
			      cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zpotrf_gpu( char uplo,  magma_int_t n, 
			      cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_zpotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
			      cuDoubleComplex *dA, magma_int_t ldda, 
			      cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              cuDoubleComplex *dA1, magma_int_t ldda1, 
                              cuDoubleComplex *dA2, magma_int_t ldda2, 
                              cuDoubleComplex *dL1, magma_int_t lddl1, 
                              cuDoubleComplex *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_ztstrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              cuDoubleComplex *hU, magma_int_t ldhu, cuDoubleComplex *dU, magma_int_t lddu, 
                              cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda, 
                              cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuDoubleComplex *hwork, magma_int_t ldhwork, cuDoubleComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_zungqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              cuDoubleComplex *da, magma_int_t ldda, 
                              cuDoubleComplex *tau, cuDoubleComplex *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_zunmqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuDoubleComplex *a,    magma_int_t lda, cuDoubleComplex *tau, 
                              cuDoubleComplex *c,    magma_int_t ldc,
                              cuDoubleComplex *work, magma_int_t lwork, 
                              cuDoubleComplex *td,   magma_int_t nb, magma_int_t *info);

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _MAGMA_Z_H_ */

