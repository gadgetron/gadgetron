/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated d
 */

#ifndef _MAGMA_D_H_
#define _MAGMA_D_H_
#define PRECISION_d

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_dgebrd( magma_int_t m, magma_int_t n, double *A, 
			  magma_int_t lda, double *d, double *e,
			  double *tauq,  double *taup, 
			  double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  double *A, magma_int_t lda, double *tau, 
			  double *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_dgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  double *A, magma_int_t lda, double *tau,
			  double *work, magma_int_t lwork,
			  double *d_T, magma_int_t *info);
magma_int_t magma_dgelqf( magma_int_t m, magma_int_t n, 
                          double *A,    magma_int_t lda,   double *tau, 
                          double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgeqlf( magma_int_t m, magma_int_t n, 
                          double *A,    magma_int_t lda,   double *tau, 
                          double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgeqrf( magma_int_t m, magma_int_t n, double *A, 
			  magma_int_t lda, double *tau, double *work, 
			  magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgetrf( magma_int_t m, magma_int_t n, double *A, 
			  magma_int_t lda, magma_int_t *ipiv, 
			  magma_int_t *info);
magma_int_t magma_dlatrd( char uplo, magma_int_t n, magma_int_t nb, double *a, 
                          magma_int_t lda, double *e, double *tau, 
			  double *w, magma_int_t ldw,
                          double *da, magma_int_t ldda, 
			  double *dw, magma_int_t lddw);
magma_int_t magma_dlahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  double *da, double *dv, double *a, 
			  magma_int_t lda, double *tau, double *t, 
			  magma_int_t ldt, double *y, magma_int_t ldy);
magma_int_t magma_dlahru( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  double *a, magma_int_t lda, 
			  double *da, double *y, 
			  double *v, double *t, 
			  double *dwork);
magma_int_t magma_dsybbd( char uplo, magma_int_t n, double *a, magma_int_t lda,
			  double *tau, double *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_dpotrf( char uplo, magma_int_t n, double *A, 
			  magma_int_t lda, magma_int_t *info);
magma_int_t magma_dsytrd( char uplo, magma_int_t n, double *A, 
			  magma_int_t lda, double *d, double *e, 
			  double *tau, double *work, magma_int_t lwork, 
			  magma_int_t *info);
magma_int_t magma_dorgqr( magma_int_t m, magma_int_t n, magma_int_t k,
			  double *a, magma_int_t lda,
			  double *tau, double *dwork,
			  magma_int_t nb, magma_int_t *info );
magma_int_t magma_dormqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          double *a, magma_int_t lda, double *tau, 
                          double *c, magma_int_t ldc, 
                          double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dormtr( char side, char uplo, char trans,
			  int m, int n,
			  double *a,    int lda,
			  double *tau,
			  double *c,    int ldc,
			  double *work, int lwork,
			  int *info);
magma_int_t magma_dorghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  double *a, magma_int_t lda,
			  double *tau,
			  double *dT, magma_int_t nb,
			  magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_dgeev( char jobvl, char jobvr, magma_int_t n,
			  double *a, magma_int_t lda,
			  double *w,
			  double *vl, magma_int_t ldvl,
			  double *vr, magma_int_t ldvr,
			  double *work, magma_int_t lwork,
			  double *rwork, magma_int_t *info);
magma_int_t magma_dgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  double *a,    magma_int_t lda, double *s, 
                          double *u,    magma_int_t ldu, 
			  double *vt,   magma_int_t ldvt,
			  double *work, magma_int_t lwork,
			  double *rwork, magma_int_t *info );
magma_int_t magma_dsyevd( char jobz, char uplo, magma_int_t n,
			  double *a, magma_int_t lda, double *w,
			  double *work, magma_int_t lwork,
		          double *rwork, magma_int_t lrwork,
		          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#else
magma_int_t  magma_dgeev( char jobvl, char jobvr, magma_int_t n,
			  double *a,    magma_int_t lda,
			  double *wr, double *wi,
			  double *vl,   magma_int_t ldvl,
			  double *vr,   magma_int_t ldvr,
			  double *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_dgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  double *a,    magma_int_t lda, double *s, 
                          double *u,    magma_int_t ldu, 
			  double *vt,   magma_int_t ldvt,
			  double *work, magma_int_t lwork,
			  magma_int_t *info );
magma_int_t magma_dsyevd( char jobz, char uplo, magma_int_t n,
			  double *a, magma_int_t lda, double *w,
			  double *work, magma_int_t lwork,
			  magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_dgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
			      double *dA,    magma_int_t ldda, 
			      double *dB,    magma_int_t lddb, 
			      double *hwork, magma_int_t lwork, 
			      magma_int_t *info);
magma_int_t magma_dgeqrf_gpu( magma_int_t m, magma_int_t n, 
			      double *dA,  magma_int_t ldda, 
			      double *tau, double *dT, 
			      magma_int_t *info);
magma_int_t magma_dgeqrf2_gpu(magma_int_t m, magma_int_t n, 
			      double *dA,  magma_int_t ldda, 
			      double *tau, magma_int_t *info);
magma_int_t magma_dgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
			      double *dA,     magma_int_t ldda, 
			      double *tau,   double *dT,
			      double *dB,    magma_int_t lddb,
			      double *hwork, magma_int_t lhwork, 
			      magma_int_t *info);
magma_int_t magma_dgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              double *dL1, magma_int_t lddl1, 
                              double *dL,  magma_int_t lddl, 
                              double *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_dgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
			      double *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dgetfl_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              double *hA, magma_int_t ldha, double *dA, magma_int_t ldda,
                              double *hL, magma_int_t ldhl, double *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              double *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_dgetrf_gpu( magma_int_t m, magma_int_t n, 
			      double *dA, magma_int_t ldda, 
			      magma_int_t *ipiv, magma_int_t *info);
magma_int_t magma_dgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
			      double *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dlabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              double *a, magma_int_t lda, double *da, magma_int_t ldda,
                              double *d, double *e, double *tauq, double *taup,  
                              double *x, magma_int_t ldx, double *dx, magma_int_t lddx, 
                              double *y, magma_int_t ldy, double *dy, magma_int_t lddy);
magma_int_t magma_dlarfb_gpu( char side, char trans, char direct, char storev, 
			      magma_int_t m, magma_int_t n, magma_int_t k,
			      double *dv, magma_int_t ldv, double *dt,    magma_int_t ldt, 
			      double *dc, magma_int_t ldc, double *dowrk, magma_int_t ldwork );
magma_int_t magma_dposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
			      double *dA, magma_int_t ldda, 
			      double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dpotrf_gpu( char uplo,  magma_int_t n, 
			      double *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_dpotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
			      double *dA, magma_int_t ldda, 
			      double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              double *dA1, magma_int_t ldda1, 
                              double *dA2, magma_int_t ldda2, 
                              double *dL1, magma_int_t lddl1, 
                              double *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_dtstrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              double *hU, magma_int_t ldhu, double *dU, magma_int_t lddu, 
                              double *hA, magma_int_t ldha, double *dA, magma_int_t ldda, 
                              double *hL, magma_int_t ldhl, double *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              double *hwork, magma_int_t ldhwork, double *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_dorgqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              double *da, magma_int_t ldda, 
                              double *tau, double *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_dormqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              double *a,    magma_int_t lda, double *tau, 
                              double *c,    magma_int_t ldc,
                              double *work, magma_int_t lwork, 
                              double *td,   magma_int_t nb, magma_int_t *info);

#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* _MAGMA_D_H_ */

