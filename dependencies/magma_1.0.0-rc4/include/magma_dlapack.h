/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated d
 */

#ifndef MAGMA_DLAPACK_H
#define MAGMA_DLAPACK_H

#define PRECISION_d
#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- LAPACK Externs used in MAGMA
*/
#if defined(ADD_)

#    define blasf77_daxpy      daxpy_
#    define blasf77_dcopy      dcopy_
#    define blasf77_ddot      ddot_ 
#    define blasf77_dgemm      dgemm_
#    define blasf77_dgemv      dgemv_
#    define blasf77_dsymm      dsymm_
#    define blasf77_dsymv      dsymv_
#    define blasf77_dsyr2k     dsyr2k_
#    define blasf77_dsyrk      dsyrk_
#    define blasf77_dscal      dscal_
#    define blasf77_dsymm      dsymm_
#    define blasf77_dsyr2k     dsyr2k_
#    define blasf77_dsyrk      dsyrk_
#    define blasf77_dswap      dswap_
#    define blasf77_dtrmm      dtrmm_
#    define blasf77_dtrmv      dtrmv_
#    define blasf77_dtrsm      dtrsm_

#    define lapackf77_dbdsqr   dbdsqr_
#    define lapackf77_dgebak   dgebak_
#    define lapackf77_dgebal   dgebal_
#    define lapackf77_dgebd2   dgebd2_
#    define lapackf77_dgebrd   dgebrd_
#    define lapackf77_dgeev    dgeev_
#    define lapackf77_dgehd2   dgehd2_
#    define lapackf77_dgehrd   dgehrd_
#    define lapackf77_dgelqf   dgelqf_
#    define lapackf77_dgels    dgels_
#    define lapackf77_dgeqlf   dgeqlf_
#    define lapackf77_dgeqrf   dgeqrf_
#    define lapackf77_dgesvd   dgesvd_
#    define lapackf77_dgetrf   dgetrf_
#    define lapackf77_dsyevd   dsyevd_
#    define lapackf77_dsytd2   dsytd2_
#    define lapackf77_dsytrd   dsytrd_
#    define lapackf77_dhseqr   dhseqr_
#    define lapackf77_dlacpy   dlacpy_
#    define lapackf77_dlacgv   dlacgv_
#    define lapackf77_dlange   dlange_
#    define lapackf77_dlansy   dlansy_
#    define lapackf77_dlansy   dlansy_
#    define lapackf77_dlarfb   dlarfb_
#    define lapackf77_dlarfg   dlarfg_
#    define lapackf77_dlarft   dlarft_
#    define lapackf77_dlarnv   dlarnv_
#    define lapackf77_dlartg   dlartg_
#    define lapackf77_dlascl   dlascl_
#    define lapackf77_dlaset   dlaset_
#    define lapackf77_dlaswp   dlaswp_
#    define lapackf77_dlatrd   dlatrd_
#    define lapackf77_dpotrf   dpotrf_
#    define lapackf77_dtrevc   dtrevc_
#    define lapackf77_dtrtri   dtrtri_
#    define lapackf77_dstedc   dstedc_
#    define lapackf77_dsymv    dsymv_
#    define lapackf77_dorg2r   dorg2r_
#    define lapackf77_dorgbr   dorgbr_
#    define lapackf77_dorghr   dorghr_
#    define lapackf77_dorglq   dorglq_
#    define lapackf77_dorgqr   dorgqr_
#    define lapackf77_dorgtr   dorgtr_
#    define lapackf77_dorm2r   dorm2r_
#    define lapackf77_dormbr   dormbr_
#    define lapackf77_dormlq   dormlq_
#    define lapackf77_dormql   dormql_
#    define lapackf77_dormqr   dormqr_
#    define lapackf77_dormtr   dormtr_

#    define lapackf77_dbdt01   dbdt01_
#    define lapackf77_dsyt21   dsyt21_
#    define lapackf77_dhst01   dhst01_
#    define lapackf77_dqrt02   dqrt02_
#    define lapackf77_dort01   dort01_
#    define lapackf77_dlarfy   dlarfy_
#    define lapackf77_dstt21   dstt21_

#elif defined(NOCHANGE)

#    define blasf77_daxpy      daxpy
#    define blasf77_dcopy      dcopy
#    define blasf77_ddot      ddot 
#    define blasf77_dgemm      dgemm
#    define blasf77_dgemv      dgemv
#    define blasf77_dsymm      dsymm
#    define blasf77_dsymv      dsymv
#    define blasf77_dsyr2k     dsyr2k
#    define blasf77_dsyrk      dsyrk
#    define blasf77_dscal      dscal
#    define blasf77_dsymm      dsymm
#    define blasf77_dsyr2k     dsyr2k
#    define blasf77_dsyrk      dsyrk
#    define blasf77_dswap      dswap
#    define blasf77_dtrmm      dtrmm
#    define blasf77_dtrmv      dtrmv
#    define blasf77_dtrsm      dtrsm

#    define lapackf77_dbdsqr   dbdsqr
#    define lapackf77_dgebak   dgebak
#    define lapackf77_dgebal   dgebal
#    define lapackf77_dgebd2   dgebd2
#    define lapackf77_dgebrd   dgebrd
#    define lapackf77_dgeev    dgeev
#    define lapackf77_dgehd2   dgehd2
#    define lapackf77_dgehrd   dgehrd
#    define lapackf77_dgelqf   dgelqf
#    define lapackf77_dgels    dgels
#    define lapackf77_dgeqlf   dgeqlf
#    define lapackf77_dgeqrf   dgeqrf
#    define lapackf77_dgesvd   dgesvd  
#    define lapackf77_dgetrf   dgetrf
#    define lapackf77_dsyevd   dsyevd
#    define lapackf77_dsytd2   dsytd2
#    define lapackf77_dsytrd   dsytrd
#    define lapackf77_dhseqr   dhseqr
#    define lapackf77_dlacpy   dlacpy
#    define lapackf77_dlacgv   dlacgv
#    define lapackf77_dlange   dlange
#    define lapackf77_dlansy   dlansy
#    define lapackf77_dlansy   dlansy
#    define lapackf77_dlarfb   dlarfb
#    define lapackf77_dlarfg   dlarfg
#    define lapackf77_dlarft   dlarft
#    define lapackf77_dlarnv   dlarnv
#    define lapackf77_dlartg   dlartg
#    define lapackf77_dlascl   dlascl
#    define lapackf77_dlaset   dlaset
#    define lapackf77_dlaswp   dlaswp
#    define lapackf77_dlatrd   dlatrd
#    define lapackf77_dpotrf   dpotrf
#    define lapackf77_dtrevc   dtrevc
#    define lapackf77_dtrtri   dtrtri
#    define lapackf77_dstedc   dstedc
#    define lapackf77_dsymv    dsymv
#    define lapackf77_dorg2r   dorg2r
#    define lapackf77_dorgbr   dorgbr
#    define lapackf77_dorghr   dorghr
#    define lapackf77_dorglq   dorglq
#    define lapackf77_dorgqr   dorgqr
#    define lapackf77_dorgtr   dorgtr
#    define lapackf77_dorm2r   dorm2r
#    define lapackf77_dormbr   dormbr
#    define lapackf77_dormlq   dormlq
#    define lapackf77_dormql   dormql
#    define lapackf77_dormqr   dormqr
#    define lapackf77_dormtr   dormtr

#    define lapackf77_dbdt01   dbdt01
#    define lapackf77_dsyt21   dsyt21
#    define lapackf77_dhst01   dhst01
#    define lapackf77_dqrt02   dqrt02
#    define lapackf77_dort01   dort01
#    define lapackf77_dlarfy   dlarfy
#    define lapackf77_dstt21   dstt21

#endif


#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        double *rwork,
#define DWORKFORZ_AND_LD double *rwork, magma_int_t *ldrwork,
#define WSPLIT           double *w
#else
#define DWORKFORZ 
#define DWORKFORZ_AND_LD
#define WSPLIT           double *wr, double *wi
#endif

  /*
   * BLAS functions (Alphabetical order)
   */
void    blasf77_daxpy( const int *, double *, double *, const int *, double *, const int *);
void    blasf77_dcopy( const int *, double *, const int *, double *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    blasf77_ddot( double *, int *, double *, int *, 
                       double *, int *);
#endif
void    blasf77_dgemm( const char *, const char *, const int *, const int *, const int *, double *, double *, const int *, double *, const int *, double *,double *, const int *);
void    blasf77_dgemv( const char *, const int  *, const int *, double *, double *, const int *, double *, const int *, double *, double *, const int *);
void    blasf77_dsymm( const char *, const char *, const int *, const int *, double *, double *, const int *, double *, const int *, double *,double *, const int *);
void    blasf77_dsymv( const char *, const int  *, double *, double *, const int *, double *, const int *, double *, double *, const int *);
void    blasf77_dsyr2k(const char *, const char *, const int *, const int *, double *, double *, const int *, double *, const int *, double *,  double *, const int *);
void    blasf77_dsyrk( const char *, const char *, const int *, const int *, double  *, double *, const int *, double *, double *, const int *);
void    blasf77_dscal( const int *, double *, double *, const int *);
void    blasf77_dsymm( const char *, const char *, const int *, const int *, double *, double *, const int *, double *, const int *, double *,double *, const int *);
void    blasf77_dsyr2k(const char *, const char *, const int *, const int *, double *, double *, const int *, double *, const int *, double *, double *, const int *);
void    blasf77_dsyrk( const char *, const char *, const int *, const int *, double *, double *, const int *, double *, double *, const int *);
void    blasf77_dswap( int *, double *, int *, double *, int *);
void    blasf77_dtrmm( const char *, const char *, const char *, const char *, const int *, const int *, double *, double *, const int *, double *,const int *);
void    blasf77_dtrmv( const char *, const char *, const char*, const int *, double*,  const int *, double *, const int*);
void    blasf77_dtrsm( const char *, const char *, const char *, const char *, const int *, const int *, double *, double *, const int *, double *,const int*);

  /*
   * Lapack functions (Alphabetical order)
   */
void    lapackf77_dbdsqr(const char *uplo, magma_int_t *n, magma_int_t *nvct, magma_int_t *nru,  magma_int_t *ncc, 
                         double *D, double *E, double *VT, magma_int_t *ldvt, double *U, magma_int_t *ldu, 
                         double *C, magma_int_t *ldc, double *work, magma_int_t *info);
void    lapackf77_dgebak(const char *job, const char *side, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         double *scale, magma_int_t *m, double *v, magma_int_t *ldv, magma_int_t *info);
void    lapackf77_dgebal(const char *job, magma_int_t *n, double *A, magma_int_t *lda, 
                         magma_int_t *ilo, magma_int_t *ihi, double *scale, magma_int_t *info);
void    lapackf77_dgebd2(magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda, double *d, double *e, double *tauq, double *taup, double *work, magma_int_t *info);
void    lapackf77_dgebrd(magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda, double *d, double *e, double *tauq, double *taup, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgeev(char *jobl, char *jobr, magma_int_t *n, double *a, magma_int_t *lda, WSPLIT, 
                        double *vl, magma_int_t *ldvl, double *vr, magma_int_t *ldvr, double *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info);
void    lapackf77_dgehd2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, double *a, magma_int_t *lda, double *tau, double *work, magma_int_t *info);
void    lapackf77_dgehrd(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, double *a, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgelqf(magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgels(const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, double *a, magma_int_t *lda, double *b, magma_int_t *ldb, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgeqlf(magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgeqrf(magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda, double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgetrf(magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info);
void    lapackf77_dgesvd(const char *jobu, const char *jobvt, magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda, double *s, double *u, magma_int_t *ldu, 
                         double *vt, magma_int_t *ldvt, double *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info );
void    lapackf77_dsyevd(const char *jobz, const char *uplo, magma_int_t *n, double *a, magma_int_t *lda, double *w, 
                         double *work, magma_int_t *lwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
void    lapackf77_dsytd2(const char *uplo, magma_int_t *n, double *a, magma_int_t *lda, double *d, double *e, double *tau, magma_int_t *info);
void    lapackf77_dsytrd(const char *uplo, magma_int_t *n, double *a, magma_int_t *lda, double *d, double *e, double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dhseqr(const char *job, const char *compz, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         double *H, magma_int_t *ldh, WSPLIT, 
                         double *Z, magma_int_t *ldz, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dlacpy(const char *uplo, magma_int_t *m, magma_int_t *n, const double *a, magma_int_t *lda, double *b, magma_int_t *ldb);
void    lapackf77_dlacgv(magma_int_t *n, double *x, magma_int_t *incx);
double  lapackf77_dlange(const char *norm, magma_int_t *m, magma_int_t *n, const double *a, magma_int_t *lda, double *work);
double  lapackf77_dlansy(const char *norm, const char *uplo, magma_int_t *n, const double *a, magma_int_t *lda, double * work);
double  lapackf77_dlansy(const char *norm, const char *uplo, magma_int_t *n, const double *a, magma_int_t *lda, double * work);
void    lapackf77_dlarfb(const char *side, const char *trans, const char *direct, const char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, const double *v, magma_int_t *ldv, const double *t, magma_int_t *ldt, double *c, magma_int_t *ldc, double *work, magma_int_t *ldwork);
void    lapackf77_dlarfg(magma_int_t *n, double *alpha, double *x, magma_int_t *incx, double *tau);
void    lapackf77_dlarft(const char *direct, const char *storev, magma_int_t *n, magma_int_t *k, double *v, magma_int_t *ldv, const double *tau, double *t, magma_int_t *ldt);
void    lapackf77_dlarnv(magma_int_t *idist, magma_int_t *iseed, magma_int_t *n, double *x);
void    lapackf77_dlartg(double *F, double *G, double *cs, double *SN, double *R);
void    lapackf77_dlascl(const char *type, magma_int_t *kl, magma_int_t *ku, double *cfrom, double *cto, 
                         magma_int_t *m, magma_int_t *n, double *A, magma_int_t *lda, magma_int_t *info);
void    lapackf77_dlaset(const char *uplo, magma_int_t *m, magma_int_t *n, double *alpha, double *beta, double *A, magma_int_t *lda);
void    lapackf77_dlaswp(magma_int_t *n, double *a, magma_int_t *lda, magma_int_t *k1, magma_int_t *k2, magma_int_t *ipiv, magma_int_t *incx);
void    lapackf77_dlatrd(const char *uplo, magma_int_t *n, magma_int_t *nb, double *a, magma_int_t *lda, double *e, double *tau, double *work, magma_int_t *ldwork);
void    lapackf77_dpotrf(const char *uplo, magma_int_t *n, double *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_dtrevc(const char *side, const char *howmny, magma_int_t *select, magma_int_t *n, 
                         double *T,  magma_int_t *ldt,  double *VL, magma_int_t *ldvl,
                         double *VR, magma_int_t *ldvr, magma_int_t *MM, magma_int_t *M, 
                         double *work, DWORKFORZ magma_int_t *info);
void    lapackf77_dtrtri(const char *uplo, const char *diag, magma_int_t *n, double *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_dstedc(const char *compz, magma_int_t *n, double *D, double *E, double *Z, magma_int_t *ldz, 
                         double *work, magma_int_t *ldwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_dsymv(const char *uplo, const magma_int_t *N, const double *alpha, const double *A, const magma_int_t *lda, const double *X, const magma_int_t *incX, const double *beta, double *Y, const magma_int_t *incY);
#endif
void    lapackf77_dorg2r(magma_int_t *m, magma_int_t *n, magma_int_t *k, double *a, magma_int_t *lda, const double *tau, double *work, magma_int_t *info);
void    lapackf77_dorgbr(const char *vect, magma_int_t *m, magma_int_t *n, magma_int_t *k, double *a, magma_int_t *lda, const double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dorghr(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, double *a, magma_int_t *lda, const double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dorglq(magma_int_t *m, magma_int_t *n, magma_int_t *k, double *a, magma_int_t *lda, const double *tau, double *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_dorgqr(magma_int_t *m, magma_int_t *n, magma_int_t *k, double *a, magma_int_t *lda, const double *tau, double *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_dorgtr(const char *uplo, magma_int_t *n, double *a, magma_int_t *lda, const double *tau, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dorm2r(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const double *a, magma_int_t *lda, const double *tau, double *c, magma_int_t *ldc, double *work, magma_int_t *info);
void    lapackf77_dormbr(const char *vect, const char *side, const char *trans, magma_int_t *M, magma_int_t *N, magma_int_t *K, double *A, magma_int_t *lda, double *Tau,
                         double *C, magma_int_t *ldc, double *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_dormlq(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const double *a, magma_int_t *lda, const double *tau, double *c, magma_int_t *ldc, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dormql(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const double *a, magma_int_t *lda, const double *tau, double *c, magma_int_t *ldc, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dormqr(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const double *a, magma_int_t *lda, const double *tau, double *c, magma_int_t *ldc, double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dormtr(const char *side, const char *uplo, const char *trans, magma_int_t *M, magma_int_t *N, double *A, magma_int_t *lda, double *Tau,
                         double *C, magma_int_t *ldc, double *work, magma_int_t *ldwork, magma_int_t *info);


  /*
   * Testing functions
   */
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_dbdt01(int *m, int *n, int *kd, double *A, int *lda, double *Q, int *ldq, double *D, double *E, double *PT, int *ldpt, double *work, double *rwork, double *resid);
void    lapackf77_dsyt21(int *itype, const char *uplo, int *n, int *kband, double *A, int *lda, double *D, double *E, double *U, int *ldu, double *V, int *ldv, double *TAU, double *work, double *rwork, double *result);
void    lapackf77_dhst01(int *n, int *ilo, int *ihi, double *A, int *lda, double *H, int *ldh, double *Q, int *ldq, double *work, int *lwork, double *rwork, double *result);
void    lapackf77_dstt21(int *n, int *kband, double *AD, double *AE, double *SD, double *SE, double *U, int *ldu, double *work, double *rwork, double *result);
void    lapackf77_dort01(const char *rowcol, int *m, int *n, double *U, int *ldu, double *work, int *lwork, double *rwork, double *resid);
#else
void    lapackf77_dbdt01(int *m, int *n, int *kd, double *A, int *lda, double *Q, int *ldq, double *D, double *E, double *PT, int *ldpt, double *work, double *resid);
void    lapackf77_dsyt21(int *itype, const char *uplo, int *n, int *kband, double *A, int *lda, double *D, double *E, double *U, int *ldu, double *V, int *ldv, double *TAU, double *work, double *result);
void    lapackf77_dhst01(int *n, int *ilo, int *ihi, double *A, int *lda, double *H, int *ldh, double *Q, int *ldq, double *work, int *lwork, double *result);
void    lapackf77_dstt21(int *n, int *kband, double *AD, double *AE, double *SD, double *SE, double *U, int *ldu, double *work, double *result);
void    lapackf77_dort01(const char *rowcol, int *m, int *n, double *U, int *ldu, double *work, int *lwork, double *resid);
#endif
void    lapackf77_dlarfy(const char *uplo, int *N, double *V, int *incv, double *tau, double *C, int *ldc, double *work);
void    lapackf77_dqrt02(int *m, int *n, int *k, double *A, double *AF, double *Q, double *R, int *lda, double *TAU, double *work, int *lwork, double *rwork, double *result);

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ 
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_d
#endif /* MAGMA ZLAPACK */
