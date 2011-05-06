/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated s
 */

#ifndef _MAGMA_S_H_
#define _MAGMA_S_H_
#define PRECISION_s

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_sgebrd( magma_int_t m, magma_int_t n, float *A, 
			  magma_int_t lda, float *d, float *e,
			  float *tauq,  float *taup, 
			  float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  float *A, magma_int_t lda, float *tau, 
			  float *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_sgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  float *A, magma_int_t lda, float *tau,
			  float *work, magma_int_t lwork,
			  float *d_T, magma_int_t *info);
magma_int_t magma_sgelqf( magma_int_t m, magma_int_t n, 
                          float *A,    magma_int_t lda,   float *tau, 
                          float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgeqlf( magma_int_t m, magma_int_t n, 
                          float *A,    magma_int_t lda,   float *tau, 
                          float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgeqrf( magma_int_t m, magma_int_t n, float *A, 
			  magma_int_t lda, float *tau, float *work, 
			  magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgetrf( magma_int_t m, magma_int_t n, float *A, 
			  magma_int_t lda, magma_int_t *ipiv, 
			  magma_int_t *info);
magma_int_t magma_slatrd( char uplo, magma_int_t n, magma_int_t nb, float *a, 
                          magma_int_t lda, float *e, float *tau, 
			  float *w, magma_int_t ldw,
                          float *da, magma_int_t ldda, 
			  float *dw, magma_int_t lddw);
magma_int_t magma_slahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  float *da, float *dv, float *a, 
			  magma_int_t lda, float *tau, float *t, 
			  magma_int_t ldt, float *y, magma_int_t ldy);
magma_int_t magma_slahru( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  float *a, magma_int_t lda, 
			  float *da, float *y, 
			  float *v, float *t, 
			  float *dwork);
magma_int_t magma_ssybbd( char uplo, magma_int_t n, float *a, magma_int_t lda,
			  float *tau, float *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_spotrf( char uplo, magma_int_t n, float *A, 
			  magma_int_t lda, magma_int_t *info);
magma_int_t magma_ssytrd( char uplo, magma_int_t n, float *A, 
			  magma_int_t lda, float *d, float *e, 
			  float *tau, float *work, magma_int_t lwork, 
			  magma_int_t *info);
magma_int_t magma_sorgqr( magma_int_t m, magma_int_t n, magma_int_t k,
			  float *a, magma_int_t lda,
			  float *tau, float *dwork,
			  magma_int_t nb, magma_int_t *info );
magma_int_t magma_sormqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          float *a, magma_int_t lda, float *tau, 
                          float *c, magma_int_t ldc, 
                          float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sormtr( char side, char uplo, char trans,
			  int m, int n,
			  float *a,    int lda,
			  float *tau,
			  float *c,    int ldc,
			  float *work, int lwork,
			  int *info);
magma_int_t magma_sorghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  float *a, magma_int_t lda,
			  float *tau,
			  float *dT, magma_int_t nb,
			  magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_sgeev( char jobvl, char jobvr, magma_int_t n,
			  float *a, magma_int_t lda,
			  float *w,
			  float *vl, magma_int_t ldvl,
			  float *vr, magma_int_t ldvr,
			  float *work, magma_int_t lwork,
			  float *rwork, magma_int_t *info);
magma_int_t magma_sgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  float *a,    magma_int_t lda, float *s, 
                          float *u,    magma_int_t ldu, 
			  float *vt,   magma_int_t ldvt,
			  float *work, magma_int_t lwork,
			  float *rwork, magma_int_t *info );
magma_int_t magma_ssyev( char jobz, char uplo, magma_int_t n,
			  float *a, magma_int_t lda, float *w,
			  float *work, magma_int_t lwork,
		          float *rwork, magma_int_t lrwork,
		          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#else
magma_int_t  magma_sgeev( char jobvl, char jobvr, magma_int_t n,
			  float *a,    magma_int_t lda,
			  float *wr, float *wi,
			  float *vl,   magma_int_t ldvl,
			  float *vr,   magma_int_t ldvr,
			  float *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_sgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  float *a,    magma_int_t lda, float *s, 
                          float *u,    magma_int_t ldu, 
			  float *vt,   magma_int_t ldvt,
			  float *work, magma_int_t lwork,
			  magma_int_t *info );
magma_int_t magma_ssyev( char jobz, char uplo, magma_int_t n,
			  float *a, magma_int_t lda, float *w,
			  float *work, magma_int_t lwork,
			  magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_sgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
			      float *dA,    magma_int_t ldda, 
			      float *dB,    magma_int_t lddb, 
			      float *hwork, magma_int_t lwork, 
			      magma_int_t *info);
magma_int_t magma_sgeqrf_gpu( magma_int_t m, magma_int_t n, 
			      float *dA,  magma_int_t ldda, 
			      float *tau, float *dT, 
			      magma_int_t *info);
magma_int_t magma_sgeqrf2_gpu(magma_int_t m, magma_int_t n, 
			      float *dA,  magma_int_t ldda, 
			      float *tau, magma_int_t *info);
magma_int_t magma_sgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
			      float *dA,     magma_int_t ldda, 
			      float *tau,   float *dT,
			      float *dB,    magma_int_t lddb,
			      float *hwork, magma_int_t lhwork, 
			      magma_int_t *info);
magma_int_t magma_sgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              float *dL1, magma_int_t lddl1, 
                              float *dL,  magma_int_t lddl, 
                              float *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_sgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
			      float *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_sgetfl_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              float *hA, magma_int_t ldha, float *dA, magma_int_t ldda,
                              float *hL, magma_int_t ldhl, float *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              float *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_sgetrf_gpu( magma_int_t m, magma_int_t n, 
			      float *dA, magma_int_t ldda, 
			      magma_int_t *ipiv, magma_int_t *info);
magma_int_t magma_sgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
			      float *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_slabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              float *a, magma_int_t lda, float *da, magma_int_t ldda,
                              float *d, float *e, float *tauq, float *taup,  
                              float *x, magma_int_t ldx, float *dx, magma_int_t lddx, 
                              float *y, magma_int_t ldy, float *dy, magma_int_t lddy);
magma_int_t magma_slarfb_gpu( char side, char trans, char direct, char storev, 
			      magma_int_t m, magma_int_t n, magma_int_t k,
			      float *dv, magma_int_t ldv, float *dt,    magma_int_t ldt, 
			      float *dc, magma_int_t ldc, float *dowrk, magma_int_t ldwork );
magma_int_t magma_sposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
			      float *dA, magma_int_t ldda, 
			      float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_spotrf_gpu( char uplo,  magma_int_t n, 
			      float *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_spotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
			      float *dA, magma_int_t ldda, 
			      float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_sssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              float *dA1, magma_int_t ldda1, 
                              float *dA2, magma_int_t ldda2, 
                              float *dL1, magma_int_t lddl1, 
                              float *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_ststrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              float *hU, magma_int_t ldhu, float *dU, magma_int_t lddu, 
                              float *hA, magma_int_t ldha, float *dA, magma_int_t ldda, 
                              float *hL, magma_int_t ldhl, float *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              float *hwork, magma_int_t ldhwork, float *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_sorgqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              float *da, magma_int_t ldda, 
                              float *tau, float *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_sormqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              float *a,    magma_int_t lda, float *tau, 
                              float *c,    magma_int_t ldc,
                              float *work, magma_int_t lwork, 
                              float *td,   magma_int_t nb, magma_int_t *info);

#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* _MAGMA_S_H_ */

