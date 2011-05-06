/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated c
 */

#ifndef MAGMA_CLAPACK_H
#define MAGMA_CLAPACK_H

#define PRECISION_c
#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- LAPACK Externs used in MAGMA
*/
#if defined(ADD_)

#    define blasf77_caxpy      caxpy_
#    define blasf77_ccopy      ccopy_
#    define blasf77_cdotc      cdotc_ 
#    define blasf77_cgemm      cgemm_
#    define blasf77_cgemv      cgemv_
#    define blasf77_chemm      chemm_
#    define blasf77_chemv      chemv_
#    define blasf77_cher2k     cher2k_
#    define blasf77_cherk      cherk_
#    define blasf77_cscal      cscal_
#    define blasf77_csymm      csymm_
#    define blasf77_csyr2k     csyr2k_
#    define blasf77_csyrk      csyrk_
#    define blasf77_cswap      cswap_
#    define blasf77_ctrmm      ctrmm_
#    define blasf77_ctrmv      ctrmv_
#    define blasf77_ctrsm      ctrsm_

#    define lapackf77_cbdsqr   cbdsqr_
#    define lapackf77_cgebak   cgebak_
#    define lapackf77_cgebal   cgebal_
#    define lapackf77_cgebd2   cgebd2_
#    define lapackf77_cgebrd   cgebrd_
#    define lapackf77_cgeev    cgeev_
#    define lapackf77_cgehd2   cgehd2_
#    define lapackf77_cgehrd   cgehrd_
#    define lapackf77_cgelqf   cgelqf_
#    define lapackf77_cgels    cgels_
#    define lapackf77_cgeqlf   cgeqlf_
#    define lapackf77_cgeqrf   cgeqrf_
#    define lapackf77_cgesvd   cgesvd_
#    define lapackf77_cgetrf   cgetrf_
#    define lapackf77_cheevd   cheevd_
#    define lapackf77_chetd2   chetd2_
#    define lapackf77_chetrd   chetrd_
#    define lapackf77_chseqr   chseqr_
#    define lapackf77_clacpy   clacpy_
#    define lapackf77_clacgv   clacgv_
#    define lapackf77_clange   clange_
#    define lapackf77_clanhe   clanhe_
#    define lapackf77_clansy   clansy_
#    define lapackf77_clarfb   clarfb_
#    define lapackf77_clarfg   clarfg_
#    define lapackf77_clarft   clarft_
#    define lapackf77_clarnv   clarnv_
#    define lapackf77_clartg   clartg_
#    define lapackf77_clascl   clascl_
#    define lapackf77_claset   claset_
#    define lapackf77_claswp   claswp_
#    define lapackf77_clatrd   clatrd_
#    define lapackf77_cpotrf   cpotrf_
#    define lapackf77_ctrevc   ctrevc_
#    define lapackf77_ctrtri   ctrtri_
#    define lapackf77_cstedc   cstedc_
#    define lapackf77_csymv    csymv_
#    define lapackf77_cung2r   cung2r_
#    define lapackf77_cungbr   cungbr_
#    define lapackf77_cunghr   cunghr_
#    define lapackf77_cunglq   cunglq_
#    define lapackf77_cungqr   cungqr_
#    define lapackf77_cungtr   cungtr_
#    define lapackf77_cunm2r   cunm2r_
#    define lapackf77_cunmbr   cunmbr_
#    define lapackf77_cunmlq   cunmlq_
#    define lapackf77_cunmql   cunmql_
#    define lapackf77_cunmqr   cunmqr_
#    define lapackf77_cunmtr   cunmtr_

#    define lapackf77_cbdt01   cbdt01_
#    define lapackf77_chet21   chet21_
#    define lapackf77_chst01   chst01_
#    define lapackf77_cqrt02   cqrt02_
#    define lapackf77_cunt01   cunt01_
#    define lapackf77_clarfy   clarfy_
#    define lapackf77_cstt21   cstt21_

#elif defined(NOCHANGE)

#    define blasf77_caxpy      caxpy
#    define blasf77_ccopy      ccopy
#    define blasf77_cdotc      cdotc 
#    define blasf77_cgemm      cgemm
#    define blasf77_cgemv      cgemv
#    define blasf77_chemm      chemm
#    define blasf77_chemv      chemv
#    define blasf77_cher2k     cher2k
#    define blasf77_cherk      cherk
#    define blasf77_cscal      cscal
#    define blasf77_csymm      csymm
#    define blasf77_csyr2k     csyr2k
#    define blasf77_csyrk      csyrk
#    define blasf77_cswap      cswap
#    define blasf77_ctrmm      ctrmm
#    define blasf77_ctrmv      ctrmv
#    define blasf77_ctrsm      ctrsm

#    define lapackf77_cbdsqr   cbdsqr
#    define lapackf77_cgebak   cgebak
#    define lapackf77_cgebal   cgebal
#    define lapackf77_cgebd2   cgebd2
#    define lapackf77_cgebrd   cgebrd
#    define lapackf77_cgeev    cgeev
#    define lapackf77_cgehd2   cgehd2
#    define lapackf77_cgehrd   cgehrd
#    define lapackf77_cgelqf   cgelqf
#    define lapackf77_cgels    cgels
#    define lapackf77_cgeqlf   cgeqlf
#    define lapackf77_cgeqrf   cgeqrf
#    define lapackf77_cgesvd   cgesvd  
#    define lapackf77_cgetrf   cgetrf
#    define lapackf77_cheevd   cheevd
#    define lapackf77_chetd2   chetd2
#    define lapackf77_chetrd   chetrd
#    define lapackf77_chseqr   chseqr
#    define lapackf77_clacpy   clacpy
#    define lapackf77_clacgv   clacgv
#    define lapackf77_clange   clange
#    define lapackf77_clanhe   clanhe
#    define lapackf77_clansy   clansy
#    define lapackf77_clarfb   clarfb
#    define lapackf77_clarfg   clarfg
#    define lapackf77_clarft   clarft
#    define lapackf77_clarnv   clarnv
#    define lapackf77_clartg   clartg
#    define lapackf77_clascl   clascl
#    define lapackf77_claset   claset
#    define lapackf77_claswp   claswp
#    define lapackf77_clatrd   clatrd
#    define lapackf77_cpotrf   cpotrf
#    define lapackf77_ctrevc   ctrevc
#    define lapackf77_ctrtri   ctrtri
#    define lapackf77_cstedc   cstedc
#    define lapackf77_csymv    csymv
#    define lapackf77_cung2r   cung2r
#    define lapackf77_cungbr   cungbr
#    define lapackf77_cunghr   cunghr
#    define lapackf77_cunglq   cunglq
#    define lapackf77_cungqr   cungqr
#    define lapackf77_cungtr   cungtr
#    define lapackf77_cunm2r   cunm2r
#    define lapackf77_cunmbr   cunmbr
#    define lapackf77_cunmlq   cunmlq
#    define lapackf77_cunmql   cunmql
#    define lapackf77_cunmqr   cunmqr
#    define lapackf77_cunmtr   cunmtr

#    define lapackf77_cbdt01   cbdt01
#    define lapackf77_chet21   chet21
#    define lapackf77_chst01   chst01
#    define lapackf77_cqrt02   cqrt02
#    define lapackf77_cunt01   cunt01
#    define lapackf77_clarfy   clarfy
#    define lapackf77_cstt21   cstt21

#endif


#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        float *rwork,
#define DWORKFORZ_AND_LD float *rwork, magma_int_t *ldrwork,
#define WSPLIT           cuFloatComplex *w
#else
#define DWORKFORZ 
#define DWORKFORZ_AND_LD
#define WSPLIT           float *wr, float *wi
#endif

  /*
   * BLAS functions (Alphabetical order)
   */
void    blasf77_caxpy( const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *);
void    blasf77_ccopy( const int *, cuFloatComplex *, const int *, cuFloatComplex *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    blasf77_cdotc( cuFloatComplex *, int *, cuFloatComplex *, int *, 
                       cuFloatComplex *, int *);
#endif
void    blasf77_cgemm( const char *, const char *, const int *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *, cuFloatComplex *,cuFloatComplex *, const int *);
void    blasf77_cgemv( const char *, const int  *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *, cuFloatComplex *, cuFloatComplex *, const int *);
void    blasf77_chemm( const char *, const char *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *, cuFloatComplex *,cuFloatComplex *, const int *);
void    blasf77_chemv( const char *, const int  *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *, cuFloatComplex *, cuFloatComplex *, const int *);
void    blasf77_cher2k(const char *, const char *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *, float *,  cuFloatComplex *, const int *);
void    blasf77_cherk( const char *, const char *, const int *, const int *, float  *, cuFloatComplex *, const int *, float *, cuFloatComplex *, const int *);
void    blasf77_cscal( const int *, cuFloatComplex *, cuFloatComplex *, const int *);
void    blasf77_csymm( const char *, const char *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *, cuFloatComplex *,cuFloatComplex *, const int *);
void    blasf77_csyr2k(const char *, const char *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, const int *, cuFloatComplex *, cuFloatComplex *, const int *);
void    blasf77_csyrk( const char *, const char *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *, cuFloatComplex *, const int *);
void    blasf77_cswap( int *, cuFloatComplex *, int *, cuFloatComplex *, int *);
void    blasf77_ctrmm( const char *, const char *, const char *, const char *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *,const int *);
void    blasf77_ctrmv( const char *, const char *, const char*, const int *, cuFloatComplex*,  const int *, cuFloatComplex *, const int*);
void    blasf77_ctrsm( const char *, const char *, const char *, const char *, const int *, const int *, cuFloatComplex *, cuFloatComplex *, const int *, cuFloatComplex *,const int*);

  /*
   * Lapack functions (Alphabetical order)
   */
void    lapackf77_cbdsqr(const char *uplo, magma_int_t *n, magma_int_t *nvct, magma_int_t *nru,  magma_int_t *ncc, 
                         float *D, float *E, cuFloatComplex *VT, magma_int_t *ldvt, cuFloatComplex *U, magma_int_t *ldu, 
                         cuFloatComplex *C, magma_int_t *ldc, float *work, magma_int_t *info);
void    lapackf77_cgebak(const char *job, const char *side, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         float *scale, magma_int_t *m, cuFloatComplex *v, magma_int_t *ldv, magma_int_t *info);
void    lapackf77_cgebal(const char *job, magma_int_t *n, cuFloatComplex *A, magma_int_t *lda, 
                         magma_int_t *ilo, magma_int_t *ihi, float *scale, magma_int_t *info);
void    lapackf77_cgebd2(magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, float *d, float *e, cuFloatComplex *tauq, cuFloatComplex *taup, cuFloatComplex *work, magma_int_t *info);
void    lapackf77_cgebrd(magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, float *d, float *e, cuFloatComplex *tauq, cuFloatComplex *taup, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgeev(char *jobl, char *jobr, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, WSPLIT, 
                        cuFloatComplex *vl, magma_int_t *ldvl, cuFloatComplex *vr, magma_int_t *ldvr, cuFloatComplex *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info);
void    lapackf77_cgehd2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *info);
void    lapackf77_cgehrd(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgelqf(magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgels(const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *b, magma_int_t *ldb, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgeqlf(magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgeqrf(magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgetrf(magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info);
void    lapackf77_cgesvd(const char *jobu, const char *jobvt, magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, float *s, cuFloatComplex *u, magma_int_t *ldu, 
                         cuFloatComplex *vt, magma_int_t *ldvt, cuFloatComplex *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info );
void    lapackf77_cheevd(const char *jobz, const char *uplo, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, float *w, 
                         cuFloatComplex *work, magma_int_t *lwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
void    lapackf77_chetd2(const char *uplo, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, float *d, float *e, cuFloatComplex *tau, magma_int_t *info);
void    lapackf77_chetrd(const char *uplo, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, float *d, float *e, cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_chseqr(const char *job, const char *compz, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         cuFloatComplex *H, magma_int_t *ldh, WSPLIT, 
                         cuFloatComplex *Z, magma_int_t *ldz, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_clacpy(const char *uplo, magma_int_t *m, magma_int_t *n, const cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *b, magma_int_t *ldb);
void    lapackf77_clacgv(magma_int_t *n, cuFloatComplex *x, magma_int_t *incx);
float  lapackf77_clange(const char *norm, magma_int_t *m, magma_int_t *n, const cuFloatComplex *a, magma_int_t *lda, float *work);
float  lapackf77_clanhe(const char *norm, const char *uplo, magma_int_t *n, const cuFloatComplex *a, magma_int_t *lda, float * work);
float  lapackf77_clansy(const char *norm, const char *uplo, magma_int_t *n, const cuFloatComplex *a, magma_int_t *lda, float * work);
void    lapackf77_clarfb(const char *side, const char *trans, const char *direct, const char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuFloatComplex *v, magma_int_t *ldv, const cuFloatComplex *t, magma_int_t *ldt, cuFloatComplex *c, magma_int_t *ldc, cuFloatComplex *work, magma_int_t *ldwork);
void    lapackf77_clarfg(magma_int_t *n, cuFloatComplex *alpha, cuFloatComplex *x, magma_int_t *incx, cuFloatComplex *tau);
void    lapackf77_clarft(const char *direct, const char *storev, magma_int_t *n, magma_int_t *k, cuFloatComplex *v, magma_int_t *ldv, const cuFloatComplex *tau, cuFloatComplex *t, magma_int_t *ldt);
void    lapackf77_clarnv(magma_int_t *idist, magma_int_t *iseed, magma_int_t *n, cuFloatComplex *x);
void    lapackf77_clartg(cuFloatComplex *F, cuFloatComplex *G, float *cs, cuFloatComplex *SN, cuFloatComplex *R);
void    lapackf77_clascl(const char *type, magma_int_t *kl, magma_int_t *ku, float *cfrom, float *cto, 
                         magma_int_t *m, magma_int_t *n, cuFloatComplex *A, magma_int_t *lda, magma_int_t *info);
void    lapackf77_claset(const char *uplo, magma_int_t *m, magma_int_t *n, cuFloatComplex *alpha, cuFloatComplex *beta, cuFloatComplex *A, magma_int_t *lda);
void    lapackf77_claswp(magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, magma_int_t *k1, magma_int_t *k2, magma_int_t *ipiv, magma_int_t *incx);
void    lapackf77_clatrd(const char *uplo, magma_int_t *n, magma_int_t *nb, cuFloatComplex *a, magma_int_t *lda, float *e, cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *ldwork);
void    lapackf77_cpotrf(const char *uplo, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_ctrevc(const char *side, const char *howmny, magma_int_t *select, magma_int_t *n, 
                         cuFloatComplex *T,  magma_int_t *ldt,  cuFloatComplex *VL, magma_int_t *ldvl,
                         cuFloatComplex *VR, magma_int_t *ldvr, magma_int_t *MM, magma_int_t *M, 
                         cuFloatComplex *work, DWORKFORZ magma_int_t *info);
void    lapackf77_ctrtri(const char *uplo, const char *diag, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_cstedc(const char *compz, magma_int_t *n, float *D, float *E, cuFloatComplex *Z, magma_int_t *ldz, 
                         cuFloatComplex *work, magma_int_t *ldwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_csymv(const char *uplo, const magma_int_t *N, const cuFloatComplex *alpha, const cuFloatComplex *A, const magma_int_t *lda, const cuFloatComplex *X, const magma_int_t *incX, const cuFloatComplex *beta, cuFloatComplex *Y, const magma_int_t *incY);
#endif
void    lapackf77_cung2r(magma_int_t *m, magma_int_t *n, magma_int_t *k, cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *info);
void    lapackf77_cungbr(const char *vect, magma_int_t *m, magma_int_t *n, magma_int_t *k, cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunghr(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunglq(magma_int_t *m, magma_int_t *n, magma_int_t *k, cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_cungqr(magma_int_t *m, magma_int_t *n, magma_int_t *k, cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_cungtr(const char *uplo, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunm2r(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc, cuFloatComplex *work, magma_int_t *info);
void    lapackf77_cunmbr(const char *vect, const char *side, const char *trans, magma_int_t *M, magma_int_t *N, magma_int_t *K, cuFloatComplex *A, magma_int_t *lda, cuFloatComplex *Tau,
                         cuFloatComplex *C, magma_int_t *ldc, cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_cunmlq(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunmql(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunmqr(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc, cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunmtr(const char *side, const char *uplo, const char *trans, magma_int_t *M, magma_int_t *N, cuFloatComplex *A, magma_int_t *lda, cuFloatComplex *Tau,
                         cuFloatComplex *C, magma_int_t *ldc, cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);


  /*
   * Testing functions
   */
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_cbdt01(int *m, int *n, int *kd, cuFloatComplex *A, int *lda, cuFloatComplex *Q, int *ldq, float *D, float *E, cuFloatComplex *PT, int *ldpt, cuFloatComplex *work, float *rwork, float *resid);
void    lapackf77_chet21(int *itype, const char *uplo, int *n, int *kband, cuFloatComplex *A, int *lda, float *D, float *E, cuFloatComplex *U, int *ldu, cuFloatComplex *V, int *ldv, cuFloatComplex *TAU, cuFloatComplex *work, float *rwork, float *result);
void    lapackf77_chst01(int *n, int *ilo, int *ihi, cuFloatComplex *A, int *lda, cuFloatComplex *H, int *ldh, cuFloatComplex *Q, int *ldq, cuFloatComplex *work, int *lwork, float *rwork, float *result);
void    lapackf77_cstt21(int *n, int *kband, float *AD, float *AE, float *SD, float *SE, cuFloatComplex *U, int *ldu, cuFloatComplex *work, float *rwork, float *result);
void    lapackf77_cunt01(const char *rowcol, int *m, int *n, cuFloatComplex *U, int *ldu, cuFloatComplex *work, int *lwork, float *rwork, float *resid);
#else
void    lapackf77_cbdt01(int *m, int *n, int *kd, cuFloatComplex *A, int *lda, cuFloatComplex *Q, int *ldq, float *D, float *E, cuFloatComplex *PT, int *ldpt, cuFloatComplex *work, float *resid);
void    lapackf77_chet21(int *itype, const char *uplo, int *n, int *kband, cuFloatComplex *A, int *lda, float *D, float *E, cuFloatComplex *U, int *ldu, cuFloatComplex *V, int *ldv, cuFloatComplex *TAU, cuFloatComplex *work, float *result);
void    lapackf77_chst01(int *n, int *ilo, int *ihi, cuFloatComplex *A, int *lda, cuFloatComplex *H, int *ldh, cuFloatComplex *Q, int *ldq, cuFloatComplex *work, int *lwork, float *result);
void    lapackf77_cstt21(int *n, int *kband, float *AD, float *AE, float *SD, float *SE, cuFloatComplex *U, int *ldu, cuFloatComplex *work, float *result);
void    lapackf77_cunt01(const char *rowcol, int *m, int *n, cuFloatComplex *U, int *ldu, cuFloatComplex *work, int *lwork, float *resid);
#endif
void    lapackf77_clarfy(const char *uplo, int *N, cuFloatComplex *V, int *incv, cuFloatComplex *tau, cuFloatComplex *C, int *ldc, cuFloatComplex *work);
void    lapackf77_cqrt02(int *m, int *n, int *k, cuFloatComplex *A, cuFloatComplex *AF, cuFloatComplex *Q, cuFloatComplex *R, int *lda, cuFloatComplex *TAU, cuFloatComplex *work, int *lwork, float *rwork, float *result);

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ 
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_c
#endif /* MAGMA ZLAPACK */
