/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated c
 */

#ifndef _MAGMA_C_H_
#define _MAGMA_C_H_
#define PRECISION_c

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_cgebrd( magma_int_t m, magma_int_t n, cuFloatComplex *A, 
			  magma_int_t lda, float *d, float *e,
			  cuFloatComplex *tauq,  cuFloatComplex *taup, 
			  cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  cuFloatComplex *A, magma_int_t lda, cuFloatComplex *tau, 
			  cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_cgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  cuFloatComplex *A, magma_int_t lda, cuFloatComplex *tau,
			  cuFloatComplex *work, magma_int_t lwork,
			  cuFloatComplex *d_T, magma_int_t *info);
magma_int_t magma_cgelqf( magma_int_t m, magma_int_t n, 
                          cuFloatComplex *A,    magma_int_t lda,   cuFloatComplex *tau, 
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgeqlf( magma_int_t m, magma_int_t n, 
                          cuFloatComplex *A,    magma_int_t lda,   cuFloatComplex *tau, 
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgeqrf( magma_int_t m, magma_int_t n, cuFloatComplex *A, 
			  magma_int_t lda, cuFloatComplex *tau, cuFloatComplex *work, 
			  magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgetrf( magma_int_t m, magma_int_t n, cuFloatComplex *A, 
			  magma_int_t lda, magma_int_t *ipiv, 
			  magma_int_t *info);
magma_int_t magma_clatrd( char uplo, magma_int_t n, magma_int_t nb, cuFloatComplex *a, 
                          magma_int_t lda, float *e, cuFloatComplex *tau, 
			  cuFloatComplex *w, magma_int_t ldw,
                          cuFloatComplex *da, magma_int_t ldda, 
			  cuFloatComplex *dw, magma_int_t lddw);
magma_int_t magma_clahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  cuFloatComplex *da, cuFloatComplex *dv, cuFloatComplex *a, 
			  magma_int_t lda, cuFloatComplex *tau, cuFloatComplex *t, 
			  magma_int_t ldt, cuFloatComplex *y, magma_int_t ldy);
magma_int_t magma_clahru( magma_int_t m, magma_int_t n, magma_int_t nb, 
			  cuFloatComplex *a, magma_int_t lda, 
			  cuFloatComplex *da, cuFloatComplex *y, 
			  cuFloatComplex *v, cuFloatComplex *t, 
			  cuFloatComplex *dwork);
magma_int_t magma_chebbd( char uplo, magma_int_t n, cuFloatComplex *a, magma_int_t lda,
			  cuFloatComplex *tau, cuFloatComplex *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_cpotrf( char uplo, magma_int_t n, cuFloatComplex *A, 
			  magma_int_t lda, magma_int_t *info);
magma_int_t magma_chetrd( char uplo, magma_int_t n, cuFloatComplex *A, 
			  magma_int_t lda, float *d, float *e, 
			  cuFloatComplex *tau, cuFloatComplex *work, magma_int_t lwork, 
			  magma_int_t *info);
magma_int_t magma_cungqr( magma_int_t m, magma_int_t n, magma_int_t k,
			  cuFloatComplex *a, magma_int_t lda,
			  cuFloatComplex *tau, cuFloatComplex *dwork,
			  magma_int_t nb, magma_int_t *info );
magma_int_t magma_cunmqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          cuFloatComplex *a, magma_int_t lda, cuFloatComplex *tau, 
                          cuFloatComplex *c, magma_int_t ldc, 
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cunmtr( char side, char uplo, char trans,
			  int m, int n,
			  cuFloatComplex *a,    int lda,
			  cuFloatComplex *tau,
			  cuFloatComplex *c,    int ldc,
			  cuFloatComplex *work, int lwork,
			  int *info);
magma_int_t magma_cunghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
			  cuFloatComplex *a, magma_int_t lda,
			  cuFloatComplex *tau,
			  cuFloatComplex *dT, magma_int_t nb,
			  magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_cgeev( char jobvl, char jobvr, magma_int_t n,
			  cuFloatComplex *a, magma_int_t lda,
			  cuFloatComplex *w,
			  cuFloatComplex *vl, magma_int_t ldvl,
			  cuFloatComplex *vr, magma_int_t ldvr,
			  cuFloatComplex *work, magma_int_t lwork,
			  float *rwork, magma_int_t *info);
magma_int_t magma_cgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  cuFloatComplex *a,    magma_int_t lda, float *s, 
                          cuFloatComplex *u,    magma_int_t ldu, 
			  cuFloatComplex *vt,   magma_int_t ldvt,
			  cuFloatComplex *work, magma_int_t lwork,
			  float *rwork, magma_int_t *info );
magma_int_t magma_cheevd( char jobz, char uplo, magma_int_t n,
			  cuFloatComplex *a, magma_int_t lda, float *w,
			  cuFloatComplex *work, magma_int_t lwork,
		          float *rwork, magma_int_t lrwork,
		          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#else
magma_int_t  magma_cgeev( char jobvl, char jobvr, magma_int_t n,
			  cuFloatComplex *a,    magma_int_t lda,
			  cuFloatComplex *wr, cuFloatComplex *wi,
			  cuFloatComplex *vl,   magma_int_t ldvl,
			  cuFloatComplex *vr,   magma_int_t ldvr,
			  cuFloatComplex *work, magma_int_t lwork,
			  magma_int_t *info);
magma_int_t magma_cgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
			  cuFloatComplex *a,    magma_int_t lda, float *s, 
                          cuFloatComplex *u,    magma_int_t ldu, 
			  cuFloatComplex *vt,   magma_int_t ldvt,
			  cuFloatComplex *work, magma_int_t lwork,
			  magma_int_t *info );
magma_int_t magma_cheevd( char jobz, char uplo, magma_int_t n,
			  cuFloatComplex *a, magma_int_t lda, float *w,
			  cuFloatComplex *work, magma_int_t lwork,
			  magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_cgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
			      cuFloatComplex *dA,    magma_int_t ldda, 
			      cuFloatComplex *dB,    magma_int_t lddb, 
			      cuFloatComplex *hwork, magma_int_t lwork, 
			      magma_int_t *info);
magma_int_t magma_cgeqrf_gpu( magma_int_t m, magma_int_t n, 
			      cuFloatComplex *dA,  magma_int_t ldda, 
			      cuFloatComplex *tau, cuFloatComplex *dT, 
			      magma_int_t *info);
magma_int_t magma_cgeqrf2_gpu(magma_int_t m, magma_int_t n, 
			      cuFloatComplex *dA,  magma_int_t ldda, 
			      cuFloatComplex *tau, magma_int_t *info);
magma_int_t magma_cgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
			      cuFloatComplex *dA,     magma_int_t ldda, 
			      cuFloatComplex *tau,   cuFloatComplex *dT,
			      cuFloatComplex *dB,    magma_int_t lddb,
			      cuFloatComplex *hwork, magma_int_t lhwork, 
			      magma_int_t *info);
magma_int_t magma_cgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              cuFloatComplex *dL1, magma_int_t lddl1, 
                              cuFloatComplex *dL,  magma_int_t lddl, 
                              cuFloatComplex *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_cgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
			      cuFloatComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_cgetfl_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              cuFloatComplex *hA, magma_int_t ldha, cuFloatComplex *dA, magma_int_t ldda,
                              cuFloatComplex *hL, magma_int_t ldhl, cuFloatComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuFloatComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_cgetrf_gpu( magma_int_t m, magma_int_t n, 
			      cuFloatComplex *dA, magma_int_t ldda, 
			      magma_int_t *ipiv, magma_int_t *info);
magma_int_t magma_cgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
			      cuFloatComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
			      cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_clabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              cuFloatComplex *a, magma_int_t lda, cuFloatComplex *da, magma_int_t ldda,
                              float *d, float *e, cuFloatComplex *tauq, cuFloatComplex *taup,  
                              cuFloatComplex *x, magma_int_t ldx, cuFloatComplex *dx, magma_int_t lddx, 
                              cuFloatComplex *y, magma_int_t ldy, cuFloatComplex *dy, magma_int_t lddy);
magma_int_t magma_clarfb_gpu( char side, char trans, char direct, char storev, 
			      magma_int_t m, magma_int_t n, magma_int_t k,
			      cuFloatComplex *dv, magma_int_t ldv, cuFloatComplex *dt,    magma_int_t ldt, 
			      cuFloatComplex *dc, magma_int_t ldc, cuFloatComplex *dowrk, magma_int_t ldwork );
magma_int_t magma_cposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
			      cuFloatComplex *dA, magma_int_t ldda, 
			      cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_cpotrf_gpu( char uplo,  magma_int_t n, 
			      cuFloatComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_cpotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
			      cuFloatComplex *dA, magma_int_t ldda, 
			      cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_cssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              cuFloatComplex *dA1, magma_int_t ldda1, 
                              cuFloatComplex *dA2, magma_int_t ldda2, 
                              cuFloatComplex *dL1, magma_int_t lddl1, 
                              cuFloatComplex *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_ctstrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              cuFloatComplex *hU, magma_int_t ldhu, cuFloatComplex *dU, magma_int_t lddu, 
                              cuFloatComplex *hA, magma_int_t ldha, cuFloatComplex *dA, magma_int_t ldda, 
                              cuFloatComplex *hL, magma_int_t ldhl, cuFloatComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuFloatComplex *hwork, magma_int_t ldhwork, cuFloatComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_cungqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              cuFloatComplex *da, magma_int_t ldda, 
                              cuFloatComplex *tau, cuFloatComplex *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_cunmqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuFloatComplex *a,    magma_int_t lda, cuFloatComplex *tau, 
                              cuFloatComplex *c,    magma_int_t ldc,
                              cuFloatComplex *work, magma_int_t lwork, 
                              cuFloatComplex *td,   magma_int_t nb, magma_int_t *info);

#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif /* _MAGMA_C_H_ */

