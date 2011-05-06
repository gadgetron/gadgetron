/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @precisions normal z -> s d c
 */

#ifndef MAGMA_ZLAPACK_H
#define MAGMA_ZLAPACK_H

#define PRECISION_z
#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- LAPACK Externs used in MAGMA
*/
#if defined(ADD_)

#    define blasf77_zaxpy      zaxpy_
#    define blasf77_zcopy      zcopy_
#    define blasf77_zdotc      zdotc_ 
#    define blasf77_zgemm      zgemm_
#    define blasf77_zgemv      zgemv_
#    define blasf77_zhemm      zhemm_
#    define blasf77_zhemv      zhemv_
#    define blasf77_zher2k     zher2k_
#    define blasf77_zherk      zherk_
#    define blasf77_zscal      zscal_
#    define blasf77_zsymm      zsymm_
#    define blasf77_zsyr2k     zsyr2k_
#    define blasf77_zsyrk      zsyrk_
#    define blasf77_zswap      zswap_
#    define blasf77_ztrmm      ztrmm_
#    define blasf77_ztrmv      ztrmv_
#    define blasf77_ztrsm      ztrsm_

#    define lapackf77_zbdsqr   zbdsqr_
#    define lapackf77_zgebak   zgebak_
#    define lapackf77_zgebal   zgebal_
#    define lapackf77_zgebd2   zgebd2_
#    define lapackf77_zgebrd   zgebrd_
#    define lapackf77_zgeev    zgeev_
#    define lapackf77_zgehd2   zgehd2_
#    define lapackf77_zgehrd   zgehrd_
#    define lapackf77_zgelqf   zgelqf_
#    define lapackf77_zgels    zgels_
#    define lapackf77_zgeqlf   zgeqlf_
#    define lapackf77_zgeqrf   zgeqrf_
#    define lapackf77_zgesvd   zgesvd_
#    define lapackf77_zgetrf   zgetrf_
#    define lapackf77_zheevd   zheevd_
#    define lapackf77_zhetd2   zhetd2_
#    define lapackf77_zhetrd   zhetrd_
#    define lapackf77_zhseqr   zhseqr_
#    define lapackf77_zlacpy   zlacpy_
#    define lapackf77_zlacgv   zlacgv_
#    define lapackf77_zlange   zlange_
#    define lapackf77_zlanhe   zlanhe_
#    define lapackf77_zlansy   zlansy_
#    define lapackf77_zlarfb   zlarfb_
#    define lapackf77_zlarfg   zlarfg_
#    define lapackf77_zlarft   zlarft_
#    define lapackf77_zlarnv   zlarnv_
#    define lapackf77_zlartg   zlartg_
#    define lapackf77_zlascl   zlascl_
#    define lapackf77_zlaset   zlaset_
#    define lapackf77_zlaswp   zlaswp_
#    define lapackf77_zlatrd   zlatrd_
#    define lapackf77_zpotrf   zpotrf_
#    define lapackf77_ztrevc   ztrevc_
#    define lapackf77_ztrtri   ztrtri_
#    define lapackf77_zstedc   zstedc_
#    define lapackf77_zsymv    zsymv_
#    define lapackf77_zung2r   zung2r_
#    define lapackf77_zungbr   zungbr_
#    define lapackf77_zunghr   zunghr_
#    define lapackf77_zunglq   zunglq_
#    define lapackf77_zungqr   zungqr_
#    define lapackf77_zungtr   zungtr_
#    define lapackf77_zunm2r   zunm2r_
#    define lapackf77_zunmbr   zunmbr_
#    define lapackf77_zunmlq   zunmlq_
#    define lapackf77_zunmql   zunmql_
#    define lapackf77_zunmqr   zunmqr_
#    define lapackf77_zunmtr   zunmtr_

#    define lapackf77_zbdt01   zbdt01_
#    define lapackf77_zhet21   zhet21_
#    define lapackf77_zhst01   zhst01_
#    define lapackf77_zqrt02   zqrt02_
#    define lapackf77_zunt01   zunt01_
#    define lapackf77_zlarfy   zlarfy_
#    define lapackf77_zstt21   zstt21_

#elif defined(NOCHANGE)

#    define blasf77_zaxpy      zaxpy
#    define blasf77_zcopy      zcopy
#    define blasf77_zdotc      zdotc 
#    define blasf77_zgemm      zgemm
#    define blasf77_zgemv      zgemv
#    define blasf77_zhemm      zhemm
#    define blasf77_zhemv      zhemv
#    define blasf77_zher2k     zher2k
#    define blasf77_zherk      zherk
#    define blasf77_zscal      zscal
#    define blasf77_zsymm      zsymm
#    define blasf77_zsyr2k     zsyr2k
#    define blasf77_zsyrk      zsyrk
#    define blasf77_zswap      zswap
#    define blasf77_ztrmm      ztrmm
#    define blasf77_ztrmv      ztrmv
#    define blasf77_ztrsm      ztrsm

#    define lapackf77_zbdsqr   zbdsqr
#    define lapackf77_zgebak   zgebak
#    define lapackf77_zgebal   zgebal
#    define lapackf77_zgebd2   zgebd2
#    define lapackf77_zgebrd   zgebrd
#    define lapackf77_zgeev    zgeev
#    define lapackf77_zgehd2   zgehd2
#    define lapackf77_zgehrd   zgehrd
#    define lapackf77_zgelqf   zgelqf
#    define lapackf77_zgels    zgels
#    define lapackf77_zgeqlf   zgeqlf
#    define lapackf77_zgeqrf   zgeqrf
#    define lapackf77_zgesvd   zgesvd  
#    define lapackf77_zgetrf   zgetrf
#    define lapackf77_zheevd   zheevd
#    define lapackf77_zhetd2   zhetd2
#    define lapackf77_zhetrd   zhetrd
#    define lapackf77_zhseqr   zhseqr
#    define lapackf77_zlacpy   zlacpy
#    define lapackf77_zlacgv   zlacgv
#    define lapackf77_zlange   zlange
#    define lapackf77_zlanhe   zlanhe
#    define lapackf77_zlansy   zlansy
#    define lapackf77_zlarfb   zlarfb
#    define lapackf77_zlarfg   zlarfg
#    define lapackf77_zlarft   zlarft
#    define lapackf77_zlarnv   zlarnv
#    define lapackf77_zlartg   zlartg
#    define lapackf77_zlascl   zlascl
#    define lapackf77_zlaset   zlaset
#    define lapackf77_zlaswp   zlaswp
#    define lapackf77_zlatrd   zlatrd
#    define lapackf77_zpotrf   zpotrf
#    define lapackf77_ztrevc   ztrevc
#    define lapackf77_ztrtri   ztrtri
#    define lapackf77_zstedc   zstedc
#    define lapackf77_zsymv    zsymv
#    define lapackf77_zung2r   zung2r
#    define lapackf77_zungbr   zungbr
#    define lapackf77_zunghr   zunghr
#    define lapackf77_zunglq   zunglq
#    define lapackf77_zungqr   zungqr
#    define lapackf77_zungtr   zungtr
#    define lapackf77_zunm2r   zunm2r
#    define lapackf77_zunmbr   zunmbr
#    define lapackf77_zunmlq   zunmlq
#    define lapackf77_zunmql   zunmql
#    define lapackf77_zunmqr   zunmqr
#    define lapackf77_zunmtr   zunmtr

#    define lapackf77_zbdt01   zbdt01
#    define lapackf77_zhet21   zhet21
#    define lapackf77_zhst01   zhst01
#    define lapackf77_zqrt02   zqrt02
#    define lapackf77_zunt01   zunt01
#    define lapackf77_zlarfy   zlarfy
#    define lapackf77_zstt21   zstt21

#endif


#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        double *rwork,
#define DWORKFORZ_AND_LD double *rwork, magma_int_t *ldrwork,
#define WSPLIT           cuDoubleComplex *w
#else
#define DWORKFORZ 
#define DWORKFORZ_AND_LD
#define WSPLIT           double *wr, double *wi
#endif

  /*
   * BLAS functions (Alphabetical order)
   */
void    blasf77_zaxpy( const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *);
void    blasf77_zcopy( const int *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    blasf77_zdotc( cuDoubleComplex *, int *, cuDoubleComplex *, int *, 
                       cuDoubleComplex *, int *);
#endif
void    blasf77_zgemm( const char *, const char *, const int *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *, cuDoubleComplex *,cuDoubleComplex *, const int *);
void    blasf77_zgemv( const char *, const int  *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *);
void    blasf77_zhemm( const char *, const char *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *, cuDoubleComplex *,cuDoubleComplex *, const int *);
void    blasf77_zhemv( const char *, const int  *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *);
void    blasf77_zher2k(const char *, const char *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *, double *,  cuDoubleComplex *, const int *);
void    blasf77_zherk( const char *, const char *, const int *, const int *, double  *, cuDoubleComplex *, const int *, double *, cuDoubleComplex *, const int *);
void    blasf77_zscal( const int *, cuDoubleComplex *, cuDoubleComplex *, const int *);
void    blasf77_zsymm( const char *, const char *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *, cuDoubleComplex *,cuDoubleComplex *, const int *);
void    blasf77_zsyr2k(const char *, const char *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *);
void    blasf77_zsyrk( const char *, const char *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *);
void    blasf77_zswap( int *, cuDoubleComplex *, int *, cuDoubleComplex *, int *);
void    blasf77_ztrmm( const char *, const char *, const char *, const char *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *,const int *);
void    blasf77_ztrmv( const char *, const char *, const char*, const int *, cuDoubleComplex*,  const int *, cuDoubleComplex *, const int*);
void    blasf77_ztrsm( const char *, const char *, const char *, const char *, const int *, const int *, cuDoubleComplex *, cuDoubleComplex *, const int *, cuDoubleComplex *,const int*);

  /*
   * Lapack functions (Alphabetical order)
   */
void    lapackf77_zbdsqr(const char *uplo, magma_int_t *n, magma_int_t *nvct, magma_int_t *nru,  magma_int_t *ncc, 
                         double *D, double *E, cuDoubleComplex *VT, magma_int_t *ldvt, cuDoubleComplex *U, magma_int_t *ldu, 
                         cuDoubleComplex *C, magma_int_t *ldc, double *work, magma_int_t *info);
void    lapackf77_zgebak(const char *job, const char *side, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         double *scale, magma_int_t *m, cuDoubleComplex *v, magma_int_t *ldv, magma_int_t *info);
void    lapackf77_zgebal(const char *job, magma_int_t *n, cuDoubleComplex *A, magma_int_t *lda, 
                         magma_int_t *ilo, magma_int_t *ihi, double *scale, magma_int_t *info);
void    lapackf77_zgebd2(magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup, cuDoubleComplex *work, magma_int_t *info);
void    lapackf77_zgebrd(magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zgeev(char *jobl, char *jobr, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, WSPLIT, 
                        cuDoubleComplex *vl, magma_int_t *ldvl, cuDoubleComplex *vr, magma_int_t *ldvr, cuDoubleComplex *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info);
void    lapackf77_zgehd2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, cuDoubleComplex *a, magma_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *info);
void    lapackf77_zgehrd(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, cuDoubleComplex *a, magma_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zgelqf(magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zgels(const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, cuDoubleComplex *a, magma_int_t *lda, cuDoubleComplex *b, magma_int_t *ldb, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zgeqlf(magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zgeqrf(magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zgetrf(magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info);
void    lapackf77_zgesvd(const char *jobu, const char *jobvt, magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, double *s, cuDoubleComplex *u, magma_int_t *ldu, 
                         cuDoubleComplex *vt, magma_int_t *ldvt, cuDoubleComplex *work, magma_int_t *lwork, DWORKFORZ magma_int_t *info );
void    lapackf77_zheevd(const char *jobz, const char *uplo, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, double *w, 
                         cuDoubleComplex *work, magma_int_t *lwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
void    lapackf77_zhetd2(const char *uplo, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, double *d, double *e, cuDoubleComplex *tau, magma_int_t *info);
void    lapackf77_zhetrd(const char *uplo, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, double *d, double *e, cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zhseqr(const char *job, const char *compz, magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         cuDoubleComplex *H, magma_int_t *ldh, WSPLIT, 
                         cuDoubleComplex *Z, magma_int_t *ldz, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zlacpy(const char *uplo, magma_int_t *m, magma_int_t *n, const cuDoubleComplex *a, magma_int_t *lda, cuDoubleComplex *b, magma_int_t *ldb);
void    lapackf77_zlacgv(magma_int_t *n, cuDoubleComplex *x, magma_int_t *incx);
double  lapackf77_zlange(const char *norm, magma_int_t *m, magma_int_t *n, const cuDoubleComplex *a, magma_int_t *lda, double *work);
double  lapackf77_zlanhe(const char *norm, const char *uplo, magma_int_t *n, const cuDoubleComplex *a, magma_int_t *lda, double * work);
double  lapackf77_zlansy(const char *norm, const char *uplo, magma_int_t *n, const cuDoubleComplex *a, magma_int_t *lda, double * work);
void    lapackf77_zlarfb(const char *side, const char *trans, const char *direct, const char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuDoubleComplex *v, magma_int_t *ldv, const cuDoubleComplex *t, magma_int_t *ldt, cuDoubleComplex *c, magma_int_t *ldc, cuDoubleComplex *work, magma_int_t *ldwork);
void    lapackf77_zlarfg(magma_int_t *n, cuDoubleComplex *alpha, cuDoubleComplex *x, magma_int_t *incx, cuDoubleComplex *tau);
void    lapackf77_zlarft(const char *direct, const char *storev, magma_int_t *n, magma_int_t *k, cuDoubleComplex *v, magma_int_t *ldv, const cuDoubleComplex *tau, cuDoubleComplex *t, magma_int_t *ldt);
void    lapackf77_zlarnv(magma_int_t *idist, magma_int_t *iseed, magma_int_t *n, cuDoubleComplex *x);
void    lapackf77_zlartg(cuDoubleComplex *F, cuDoubleComplex *G, double *cs, cuDoubleComplex *SN, cuDoubleComplex *R);
void    lapackf77_zlascl(const char *type, magma_int_t *kl, magma_int_t *ku, double *cfrom, double *cto, 
                         magma_int_t *m, magma_int_t *n, cuDoubleComplex *A, magma_int_t *lda, magma_int_t *info);
void    lapackf77_zlaset(const char *uplo, magma_int_t *m, magma_int_t *n, cuDoubleComplex *alpha, cuDoubleComplex *beta, cuDoubleComplex *A, magma_int_t *lda);
void    lapackf77_zlaswp(magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, magma_int_t *k1, magma_int_t *k2, magma_int_t *ipiv, magma_int_t *incx);
void    lapackf77_zlatrd(const char *uplo, magma_int_t *n, magma_int_t *nb, cuDoubleComplex *a, magma_int_t *lda, double *e, cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *ldwork);
void    lapackf77_zpotrf(const char *uplo, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_ztrevc(const char *side, const char *howmny, magma_int_t *select, magma_int_t *n, 
                         cuDoubleComplex *T,  magma_int_t *ldt,  cuDoubleComplex *VL, magma_int_t *ldvl,
                         cuDoubleComplex *VR, magma_int_t *ldvr, magma_int_t *MM, magma_int_t *M, 
                         cuDoubleComplex *work, DWORKFORZ magma_int_t *info);
void    lapackf77_ztrtri(const char *uplo, const char *diag, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_zstedc(const char *compz, magma_int_t *n, double *D, double *E, cuDoubleComplex *Z, magma_int_t *ldz, 
                         cuDoubleComplex *work, magma_int_t *ldwork, DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_zsymv(const char *uplo, const magma_int_t *N, const cuDoubleComplex *alpha, const cuDoubleComplex *A, const magma_int_t *lda, const cuDoubleComplex *X, const magma_int_t *incX, const cuDoubleComplex *beta, cuDoubleComplex *Y, const magma_int_t *incY);
#endif
void    lapackf77_zung2r(magma_int_t *m, magma_int_t *n, magma_int_t *k, cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *info);
void    lapackf77_zungbr(const char *vect, magma_int_t *m, magma_int_t *n, magma_int_t *k, cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zunghr(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zunglq(magma_int_t *m, magma_int_t *n, magma_int_t *k, cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_zungqr(magma_int_t *m, magma_int_t *n, magma_int_t *k, cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_zungtr(const char *uplo, magma_int_t *n, cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zunm2r(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *c, magma_int_t *ldc, cuDoubleComplex *work, magma_int_t *info);
void    lapackf77_zunmbr(const char *vect, const char *side, const char *trans, magma_int_t *M, magma_int_t *N, magma_int_t *K, cuDoubleComplex *A, magma_int_t *lda, cuDoubleComplex *Tau,
                         cuDoubleComplex *C, magma_int_t *ldc, cuDoubleComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_zunmlq(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *c, magma_int_t *ldc, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zunmql(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *c, magma_int_t *ldc, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zunmqr(const char *side, const char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *k, const cuDoubleComplex *a, magma_int_t *lda, const cuDoubleComplex *tau, cuDoubleComplex *c, magma_int_t *ldc, cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_zunmtr(const char *side, const char *uplo, const char *trans, magma_int_t *M, magma_int_t *N, cuDoubleComplex *A, magma_int_t *lda, cuDoubleComplex *Tau,
                         cuDoubleComplex *C, magma_int_t *ldc, cuDoubleComplex *work, magma_int_t *ldwork, magma_int_t *info);


  /*
   * Testing functions
   */
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_zbdt01(int *m, int *n, int *kd, cuDoubleComplex *A, int *lda, cuDoubleComplex *Q, int *ldq, double *D, double *E, cuDoubleComplex *PT, int *ldpt, cuDoubleComplex *work, double *rwork, double *resid);
void    lapackf77_zhet21(int *itype, const char *uplo, int *n, int *kband, cuDoubleComplex *A, int *lda, double *D, double *E, cuDoubleComplex *U, int *ldu, cuDoubleComplex *V, int *ldv, cuDoubleComplex *TAU, cuDoubleComplex *work, double *rwork, double *result);
void    lapackf77_zhst01(int *n, int *ilo, int *ihi, cuDoubleComplex *A, int *lda, cuDoubleComplex *H, int *ldh, cuDoubleComplex *Q, int *ldq, cuDoubleComplex *work, int *lwork, double *rwork, double *result);
void    lapackf77_zstt21(int *n, int *kband, double *AD, double *AE, double *SD, double *SE, cuDoubleComplex *U, int *ldu, cuDoubleComplex *work, double *rwork, double *result);
void    lapackf77_zunt01(const char *rowcol, int *m, int *n, cuDoubleComplex *U, int *ldu, cuDoubleComplex *work, int *lwork, double *rwork, double *resid);
#else
void    lapackf77_zbdt01(int *m, int *n, int *kd, cuDoubleComplex *A, int *lda, cuDoubleComplex *Q, int *ldq, double *D, double *E, cuDoubleComplex *PT, int *ldpt, cuDoubleComplex *work, double *resid);
void    lapackf77_zhet21(int *itype, const char *uplo, int *n, int *kband, cuDoubleComplex *A, int *lda, double *D, double *E, cuDoubleComplex *U, int *ldu, cuDoubleComplex *V, int *ldv, cuDoubleComplex *TAU, cuDoubleComplex *work, double *result);
void    lapackf77_zhst01(int *n, int *ilo, int *ihi, cuDoubleComplex *A, int *lda, cuDoubleComplex *H, int *ldh, cuDoubleComplex *Q, int *ldq, cuDoubleComplex *work, int *lwork, double *result);
void    lapackf77_zstt21(int *n, int *kband, double *AD, double *AE, double *SD, double *SE, cuDoubleComplex *U, int *ldu, cuDoubleComplex *work, double *result);
void    lapackf77_zunt01(const char *rowcol, int *m, int *n, cuDoubleComplex *U, int *ldu, cuDoubleComplex *work, int *lwork, double *resid);
#endif
void    lapackf77_zlarfy(const char *uplo, int *N, cuDoubleComplex *V, int *incv, cuDoubleComplex *tau, cuDoubleComplex *C, int *ldc, cuDoubleComplex *work);
void    lapackf77_zqrt02(int *m, int *n, int *k, cuDoubleComplex *A, cuDoubleComplex *AF, cuDoubleComplex *Q, cuDoubleComplex *R, int *lda, cuDoubleComplex *TAU, cuDoubleComplex *work, int *lwork, double *rwork, double *result);

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ 
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_z
#endif /* MAGMA ZLAPACK */
