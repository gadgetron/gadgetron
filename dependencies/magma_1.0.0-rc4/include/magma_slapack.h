/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated s
 */

#ifndef MAGMA_SLAPACK_H
#define MAGMA_SLAPACK_H

#define PRECISION_s
#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- LAPACK Externs used in MAGMA
*/
#if defined(ADD_)

#    define blasf77_saxpy      saxpy_
#    define blasf77_scopy      scopy_
#    define blasf77_sdot      sdot_ 
#    define blasf77_sgemm      sgemm_
#    define blasf77_sgemv      sgemv_
#    define blasf77_ssymm      ssymm_
#    define blasf77_ssymv      ssymv_
#    define blasf77_ssyr2k     ssyr2k_
#    define blasf77_ssyrk      ssyrk_
#    define blasf77_sscal      sscal_
#    define blasf77_ssymm      ssymm_
#    define blasf77_ssyr2k     ssyr2k_
#    define blasf77_ssyrk      ssyrk_
#    define blasf77_sswap      sswap_
#    define blasf77_strmm      strmm_
#    define blasf77_strmv      strmv_
#    define blasf77_strsm      strsm_

#    define lapackf77_sbdsqr   sbdsqr_
#    define lapackf77_sgebak   sgebak_
#    define lapackf77_sgebal   sgebal_
#    define lapackf77_sgebd2   sgebd2_
#    define lapackf77_sgebrd   sgebrd_
#    define lapackf77_sgeev    sgeev_
#    define lapackf77_sgehd2   sgehd2_
#    define lapackf77_sgehrd   sgehrd_
#    define lapackf77_sgelqf   sgelqf_
#    define lapackf77_sgels    sgels_
#    define lapackf77_sgeqlf   sgeqlf_
#    define lapackf77_sgeqrf   sgeqrf_
#    define lapackf77_sgesvd   sgesvd_
#    define lapackf77_sgetrf   sgetrf_
#    define lapackf77_ssyev   ssyev_
#    define lapackf77_ssytd2   ssytd2_
#    define lapackf77_ssytrd   ssytrd_
#    define lapackf77_shseqr   shseqr_
#    define lapackf77_slacpy   slacpy_
#    define lapackf77_slacgv   slacgv_
#    define lapackf77_slange   slange_
#    define lapackf77_slansy   slansy_
#    define lapackf77_slansy   slansy_
#    define lapackf77_slarfb   slarfb_
#    define lapackf77_slarfg   slarfg_
#    define lapackf77_slarft   slarft_
#    define lapackf77_slarnv   slarnv_
#    define lapackf77_slartg   slartg_
#    define lapackf77_slascl   slascl_
#    define lapackf77_slaset   slaset_
#    define lapackf77_slaswp   slaswp_
#    define lapackf77_slatrd   slatrd_
#    define lapackf77_spotrf   spotrf_
#    define lapackf77_strevc   strevc_
#    define lapackf77_strtri   strtri_
#    define lapackf77_sstedc   sstedc_
#    define lapackf77_ssymv    ssymv_
#    define lapackf77_sorg2r   sorg2r_
#    define lapackf77_sorgbr   sorgbr_
#    define lapackf77_sorghr   sorghr_
#    define lapackf77_sorglq   sorglq_
#    define lapackf77_sorgqr   sorgqr_
#    define lapackf77_sorgtr   sorgtr_
#    define lapackf77_sorm2r   sorm2r_
#    define lapackf77_sormbr   sormbr_
#    define lapackf77_sormlq   sormlq_
#    define lapackf77_sormql   sormql_
#    define lapackf77_sormqr   sormqr_
#    define lapackf77_sormtr   sormtr_

#    define lapackf77_sbdt01   sbdt01_
#    define lapackf77_ssyt21   ssyt21_
#    define lapackf77_shst01   shst01_
#    define lapackf77_sqrt02   sqrt02_
#    define lapackf77_sort01   sort01_
#    define lapackf77_slarfy   slarfy_
#    define lapackf77_sstt21   sstt21_

#elif defined(NOCHANGE)

#    define blasf77_saxpy      saxpy
#    define blasf77_scopy      scopy
#    define blasf77_sdot      sdot 
#    define blasf77_sgemm      sgemm
#    define blasf77_sgemv      sgemv
#    define blasf77_ssymm      ssymm
#    define blasf77_ssymv      ssymv
#    define blasf77_ssyr2k     ssyr2k
#    define blasf77_ssyrk      ssyrk
#    define blasf77_sscal      sscal
#    define blasf77_ssymm      ssymm
#    define blasf77_ssyr2k     ssyr2k
#    define blasf77_ssyrk      ssyrk
#    define blasf77_sswap      sswap
#    define blasf77_strmm      strmm
#    define blasf77_strmv      strmv
#    define blasf77_strsm      strsm

#    define lapackf77_sbdsqr   sbdsqr
#    define lapackf77_sgebak   sgebak
#    define lapackf77_sgebal   sgebal
#    define lapackf77_sgebd2   sgebd2
#    define lapackf77_sgebrd   sgebrd
#    define lapackf77_sgeev    sgeev
#    define lapackf77_sgehd2   sgehd2
#    define lapackf77_sgehrd   sgehrd
#    define lapackf77_sgelqf   sgelqf
#    define lapackf77_sgels    sgels
#    define lapackf77_sgeqlf   sgeqlf
#    define lapackf77_sgeqrf   sgeqrf
#    define lapackf77_sgesvd   sgesvd  
#    define lapackf77_sgetrf   sgetrf
#    define lapackf77_ssyev   ssyev
#    define lapackf77_ssytd2   ssytd2
#    define lapackf77_ssytrd   ssytrd
#    define lapackf77_shseqr   shseqr
#    define lapackf77_slacpy   slacpy
#    define lapackf77_slacgv   slacgv
#    define lapackf77_slange   slange
#    define lapackf77_slansy   slansy
#    define lapackf77_slansy   slansy
#    define lapackf77_slarfb   slarfb
#    define lapackf77_slarfg   slarfg
#    define lapackf77_slarft   slarft
#    define lapackf77_slarnv   slarnv
#    define lapackf77_slartg   slartg
#    define lapackf77_slascl   slascl
#    define lapackf77_slaset   slaset
#    define lapackf77_slaswp   slaswp
#    define lapackf77_slatrd   slatrd
#    define lapackf77_spotrf   spotrf
#    define lapackf77_strevc   strevc
#    define lapackf77_strtri   strtri
#    define lapackf77_sstedc   sstedc
#    define lapackf77_ssymv    ssymv
#    define lapackf77_sorg2r   sorg2r
#    define lapackf77_sorgbr   sorgbr
#    define lapackf77_sorghr   sorghr
#    define lapackf77_sorglq   sorglq
#    define lapackf77_sorgqr   sorgqr
#    define lapackf77_sorgtr   sorgtr
#    define lapackf77_sorm2r   sorm2r
#    define lapackf77_sormbr   sormbr
#    define lapackf77_sormlq   sormlq
#    define lapackf77_sormql   sormql
#    define lapackf77_sormqr   sormqr
#    define lapackf77_sormtr   sormtr

#    define lapackf77_sbdt01   sbdt01
#    define lapackf77_ssyt21   ssyt21
#    define lapackf77_shst01   shst01
#    define lapackf77_sqrt02   sqrt02
#    define lapackf77_sort01   sort01
#    define lapackf77_slarfy   slarfy
#    define lapackf77_sstt21   sstt21

#endif


#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        float *rwork,
#define DWORKFORZ_AND_LD float *rwork, magma_int_t *ldrwork,
#define WSPLIT           float *w
#else
#define DWORKFORZ 
#define DWORKFORZ_AND_LD
#define WSPLIT           float *wr, float *wi
#endif

  /*
   * BLAS functions (Alphabetical order)
   */
void    blasf77_saxpy( const int *, float *, float *, const int *, float *, const int *);
void    blasf77_scopy( const int *, float *, const int *, float *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    blasf77_sdot( float *, int *, float *, int *, 
                       float *, int *);
#endif
void    blasf77_sgemm( const char *, const char *, const int *, const int *, const int *, float *, float *, const int *, float *, const int *, float *,float *, const int *);
void    blasf77_sgemv( const char *, const int  *, const int *, float *, float *, const int *, float *, const int *, float *, float *, const int *);
void    blasf77_ssymm( const char *, const char *, const int *, const int *, float *, float *, const int *, float *, const int *, float *,float *, const int *);
void    blasf77_ssymv( const char *, const int  *, float *, float *, const int *, float *, const int *, float *, float *, const int *);
void    blasf77_ssyr2k(const char *, const char *, const int *, const int *, float *, float *, const int *, float *, const int *, float *,  float *, const int *);
void    blasf77_ssyrk( const char *, const char *, const int *, const int *, float  *, float *, const int *, float *, float *, const int *);
void    blasf77_sscal( const int *, float *, float *, const int *);
void    blasf77_ssymm( const char *, const char *, const int *, const int *, float *, float *, const int *, float *, const int *, float *,float *, const int *);
void    blasf77_ssyr2k(const char *, const char *, const int *, const int *, float *, float *, const int *, float *, const int *, float *, float *, const int *);
void    blasf77_ssyrk( const char *, const char *, const int *, const int *, float *, float *, const int *, float *, float *, const int *);
void    blasf77_sswap( int *, float *, int *, float *, int *);
void    blasf77_strmm( const char *, const char *, const char *, const char *, const int *, const int *, float *, float *, const int *, float *,const int *);
void    blasf77_strmv( const char *, const char *, const char*, const int *, float*,  const int *, float *, const int*);
void    blasf77_strsm( const char *, const char *, const char *, const char *, const int *, const int *, float *, float *, const int *, float *,const int*);

  /*
   * Lapack functions (Alphabetical order)
   */
void    lapackf77_sbdsqr(const char *uplo, magma_int_t *n, magma_int_t *nvct, magma_int_t *nru,  magma_int_t *ncc, 
                         float *D, float *E, float *VT, magma_int_t *ldvt, float *U, magma_int_t *ldu, 
                         float *C, magma_int_t *ldc, float *work, magma_int_t *info);
void    lapackf77_sgebak(const char *job, const char *side, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         float *scale, magma_int_t *m, float *v, magma_int_t *ldv, magma_int_t *info);
void    lapackf77_sgebal(const char *job, magma_int_t *n, float *A, magma_int_t *lda, 
                         magma_int_t *ilo, magma_int_t *ihi, float *scale, magma_int_t *info);
void    lapackf77_sgebd2(magma_int_t *m, magma_int_t *n, float *a, magma_int_t *lda, float *d, float *e, float *tauq, float *taup, float *work, magma_int_t *info);
void    lapackf77_sgebrd(magma_int_t *m, magma_int_t *n, float *a, magma_int_t *lda, float *d, float *e, float *tauq, float *taup, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sgeev(char *jobl, char *jobr, magma_int_t *n, float *a, magma_int_t *lda, WSPLIT, 
                        float *vl, magma_int_t *ldvl, float *vr, magma_int_t *ldvr, float *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info);
void    lapackf77_sgehd2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, float *a, magma_int_t *lda, float *tau, float *work, magma_int_t *info);
void    lapackf77_sgehrd(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, float *a, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sgelqf(magma_int_t *m, magma_int_t *n, float *a, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sgels(const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, float *a, magma_int_t *lda, float *b, magma_int_t *ldb, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sgeqlf(magma_int_t *m, magma_int_t *n, float *a, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sgeqrf(magma_int_t *m, magma_int_t *n, float *a, magma_int_t *lda, float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sgetrf(magma_int_t *m, magma_int_t *n, float *a, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info);
void    lapackf77_sgesvd(const char *jobu, const char *jobvt, magma_int_t *m, magma_int_t *n, float *a, magma_int_t *lda, float *s, float *u, magma_int_t *ldu, 
                         float *vt, magma_int_t *ldvt, float *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info );
void    lapackf77_ssyev(const char *jobz, const char *uplo, magma_int_t *n, float *a, magma_int_t *lda, float *w, 
                         float *work, magma_int_t *lwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
void    lapackf77_ssytd2(const char *uplo, magma_int_t *n, float *a, magma_int_t *lda, float *d, float *e, float *tau, magma_int_t *info);
void    lapackf77_ssytrd(const char *uplo, magma_int_t *n, float *a, magma_int_t *lda, float *d, float *e, float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_shseqr(const char *job, const char *compz, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         float *H, magma_int_t *ldh, WSPLIT, 
                         float *Z, magma_int_t *ldz, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_slacpy(const char *uplo, magma_int_t *m, magma_int_t *n, const float *a, magma_int_t *lda, float *b, magma_int_t *ldb);
void    lapackf77_slacgv(magma_int_t *n, float *x, magma_int_t *incx);
float  lapackf77_slange(const char *norm, magma_int_t *m, magma_int_t *n, const float *a, magma_int_t *lda, float *work);
float  lapackf77_slansy(const char *norm, const char *uplo, magma_int_t *n, const float *a, magma_int_t *lda, float * work);
float  lapackf77_slansy(const char *norm, const char *uplo, magma_int_t *n, const float *a, magma_int_t *lda, float * work);
void    lapackf77_slarfb(const char *side, const char *trans, const char *direct, const char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, const float *v, magma_int_t *ldv, const float *t, magma_int_t *ldt, float *c, magma_int_t *ldc, float *work, magma_int_t *ldwork);
void    lapackf77_slarfg(magma_int_t *n, float *alpha, float *x, magma_int_t *incx, float *tau);
void    lapackf77_slarft(const char *direct, const char *storev, magma_int_t *n, magma_int_t *k, float *v, magma_int_t *ldv, const float *tau, float *t, magma_int_t *ldt);
void    lapackf77_slarnv(magma_int_t *idist, magma_int_t *iseed, magma_int_t *n, float *x);
void    lapackf77_slartg(float *F, float *G, float *cs, float *SN, float *R);
void    lapackf77_slascl(const char *type, magma_int_t *kl, magma_int_t *ku, float *cfrom, float *cto, 
                         magma_int_t *m, magma_int_t *n, float *A, magma_int_t *lda, magma_int_t *info);
void    lapackf77_slaset(const char *uplo, magma_int_t *m, magma_int_t *n, float *alpha, float *beta, float *A, magma_int_t *lda);
void    lapackf77_slaswp(magma_int_t *n, float *a, magma_int_t *lda, magma_int_t *k1, magma_int_t *k2, magma_int_t *ipiv, magma_int_t *incx);
void    lapackf77_slatrd(const char *uplo, magma_int_t *n, magma_int_t *nb, float *a, magma_int_t *lda, float *e, float *tau, float *work, magma_int_t *ldwork);
void    lapackf77_spotrf(const char *uplo, magma_int_t *n, float *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_strevc(const char *side, const char *howmny, magma_int_t *select, magma_int_t *n, 
                         float *T,  magma_int_t *ldt,  float *VL, magma_int_t *ldvl,
                         float *VR, magma_int_t *ldvr, magma_int_t *MM, magma_int_t *M, 
                         float *work, DWORKFORZ magma_int_t *info);
void    lapackf77_strtri(const char *uplo, const char *diag, magma_int_t *n, float *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_sstedc(const char *compz, magma_int_t *n, float *D, float *E, float *Z, magma_int_t *ldz, 
                         float *work, magma_int_t *ldwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_ssymv(const char *uplo, const magma_int_t *N, const float *alpha, const float *A, const magma_int_t *lda, const float *X, const magma_int_t *incX, const float *beta, float *Y, const magma_int_t *incY);
#endif
void    lapackf77_sorg2r(magma_int_t *m, magma_int_t *n, magma_int_t *k, float *a, magma_int_t *lda, const float *tau, float *work, magma_int_t *info);
void    lapackf77_sorgbr(const char *vect, magma_int_t *m, magma_int_t *n, magma_int_t *k, float *a, magma_int_t *lda, const float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sorghr(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, float *a, magma_int_t *lda, const float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sorglq(magma_int_t *m, magma_int_t *n, magma_int_t *k, float *a, magma_int_t *lda, const float *tau, float *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_sorgqr(magma_int_t *m, magma_int_t *n, magma_int_t *k, float *a, magma_int_t *lda, const float *tau, float *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_sorgtr(const char *uplo, magma_int_t *n, float *a, magma_int_t *lda, const float *tau, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sorm2r(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const float *a, magma_int_t *lda, const float *tau, float *c, magma_int_t *ldc, float *work, magma_int_t *info);
void    lapackf77_sormbr(const char *vect, const char *side, const char *trans, magma_int_t *M, magma_int_t *N, magma_int_t *K, float *A, magma_int_t *lda, float *Tau,
                         float *C, magma_int_t *ldc, float *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_sormlq(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const float *a, magma_int_t *lda, const float *tau, float *c, magma_int_t *ldc, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sormql(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const float *a, magma_int_t *lda, const float *tau, float *c, magma_int_t *ldc, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sormqr(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const float *a, magma_int_t *lda, const float *tau, float *c, magma_int_t *ldc, float *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_sormtr(const char *side, const char *uplo, const char *trans, magma_int_t *M, magma_int_t *N, float *A, magma_int_t *lda, float *Tau,
                         float *C, magma_int_t *ldc, float *work, magma_int_t *ldwork, magma_int_t *info);


  /*
   * Testing functions
   */
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_sbdt01(int *m, int *n, int *kd, float *A, int *lda, float *Q, int *ldq, float *D, float *E, float *PT, int *ldpt, float *work, float *rwork, float *resid);
void    lapackf77_ssyt21(int *itype, const char *uplo, int *n, int *kband, float *A, int *lda, float *D, float *E, float *U, int *ldu, float *V, int *ldv, float *TAU, float *work, float *rwork, float *result);
void    lapackf77_shst01(int *n, int *ilo, int *ihi, float *A, int *lda, float *H, int *ldh, float *Q, int *ldq, float *work, int *lwork, float *rwork, float *result);
void    lapackf77_sstt21(int *n, int *kband, float *AD, float *AE, float *SD, float *SE, float *U, int *ldu, float *work, float *rwork, float *result);
void    lapackf77_sort01(const char *rowcol, int *m, int *n, float *U, int *ldu, float *work, int *lwork, float *rwork, float *resid);
#else
void    lapackf77_sbdt01(int *m, int *n, int *kd, float *A, int *lda, float *Q, int *ldq, float *D, float *E, float *PT, int *ldpt, float *work, float *resid);
void    lapackf77_ssyt21(int *itype, const char *uplo, int *n, int *kband, float *A, int *lda, float *D, float *E, float *U, int *ldu, float *V, int *ldv, float *TAU, float *work, float *result);
void    lapackf77_shst01(int *n, int *ilo, int *ihi, float *A, int *lda, float *H, int *ldh, float *Q, int *ldq, float *work, int *lwork, float *result);
void    lapackf77_sstt21(int *n, int *kband, float *AD, float *AE, float *SD, float *SE, float *U, int *ldu, float *work, float *result);
void    lapackf77_sort01(const char *rowcol, int *m, int *n, float *U, int *ldu, float *work, int *lwork, float *resid);
#endif
void    lapackf77_slarfy(const char *uplo, int *N, float *V, int *incv, float *tau, float *C, int *ldc, float *work);
void    lapackf77_sqrt02(int *m, int *n, int *k, float *A, float *AF, float *Q, float *R, int *lda, float *TAU, float *work, int *lwork, float *rwork, float *result);

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ 
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_s
#endif /* MAGMA ZLAPACK */
