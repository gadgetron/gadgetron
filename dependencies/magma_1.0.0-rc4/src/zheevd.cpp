/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> c

*/
#include "common_magma.h"

extern "C" magma_int_t 
magma_zheevd(char jobz, char uplo, 
	     magma_int_t n, 
	     cuDoubleComplex *a, magma_int_t lda, 
	     double *w, 
	     cuDoubleComplex *work, magma_int_t lwork,
	     double *rwork, magma_int_t lrwork,
	     magma_int_t *iwork, magma_int_t liwork,
	     magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======
    ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a   
    complex Hermitian matrix A.  If eigenvectors are desired, it uses a   
    divide and conquer algorithm.   

    The divide and conquer algorithm makes very mild assumptions about   
    floating point arithmetic. It will work on machines with a guard   
    digit in add/subtract, or on those binary machines without guard   
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or   
    Cray-2. It could conceivably fail on hexadecimal or decimal machines   
    without guard digits, but we know of none.   

    Arguments   
    =========   
    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX_16 array, dimension (LDA, N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')   
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) COMPLEX_16 array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.   
            If N <= 1,                LWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.   
            If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal sizes of the WORK, RWORK and   
            IWORK arrays, returns these values as the first entries of   
            the WORK, RWORK and IWORK arrays, and no error message   
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.    

    RWORK   (workspace/output) DOUBLE PRECISION array,   
                                           dimension (LRWORK)   
            On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.   

    LRWORK  (input) INTEGER   
            The dimension of the array RWORK.   
            If N <= 1,                LRWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LRWORK must be at least N.   
            If JOBZ  = 'V' and N > 1, LRWORK must be at least   
                           1 + 5*N + 2*N**2.   

            If LRWORK = -1, then a workspace query is assumed; the   
            routine only calculates the optimal sizes of the WORK, RWORK   
            and IWORK arrays, returns these values as the first entries   
            of the WORK, RWORK and IWORK arrays, and no error message   
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.   

    IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))   
            On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.   

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If N <= 1,                LIWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.   

            If LIWORK = -1, then a workspace query is assumed; the   
            routine only calculates the optimal sizes of the WORK, RWORK   
            and IWORK arrays, returns these values as the first entries   
            of the WORK, RWORK and IWORK arrays, and no error message   
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed   
                  to converge; i off-diagonal elements of an intermediate   
                  tridiagonal form did not converge to zero;   
                  if INFO = i and JOBZ = 'V', then the algorithm failed   
                  to compute an eigenvalue while working on the submatrix   
                  lying in rows and columns INFO/(N+1) through   
                  mod(INFO,N+1).   

    Further Details   
    ===============   
    Based on contributions by   
       Jeff Rutter, Computer Science Division, University of California   
       at Berkeley, USA   

    Modified description of INFO. Sven, 16 Feb 05.   
    =====================================================================   */


    char uplo_[2] = {uplo, 0};
    char jobz_[2] = {jobz, 0};
    static magma_int_t c__1 = 1;
    static magma_int_t c_n1 = -1;
    static magma_int_t c__0 = 0;
    static double c_b18 = 1.;
    
    magma_int_t a_dim1, a_offset;
    double d__1;

    static double eps;
    static magma_int_t inde;
    static double anrm;
    static magma_int_t imax;
    static double rmin, rmax;
    static double sigma;
    static magma_int_t iinfo, lwmin;
    static magma_int_t lower;
    static magma_int_t llrwk;
    static magma_int_t wantz;
    static magma_int_t indwk2, llwrk2;
    static magma_int_t iscale;
    static double safmin;
    static double bignum;
    static magma_int_t indtau;
    static magma_int_t indrwk, indwrk, liwmin;
    static magma_int_t lrwmin, llwork;
    static double smlnum;
    static magma_int_t lquery;

    /* Function Body */
    wantz = lapackf77_lsame(jobz_, MagmaVectorsStr);
    lower = lapackf77_lsame(uplo_, MagmaLowerStr);
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    *info = 0;
    if (! (wantz || lapackf77_lsame(jobz_, MagmaNoVectorsStr))) {
	*info = -1;
    } else if (! (lower || lapackf77_lsame(uplo_, MagmaUpperStr))) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max(1,n)) {
	*info = -5;
    }

    lapackf77_zheevd(jobz_, uplo_, &n, 
                     a, &lda, w, work, &c_n1, 
                     rwork, &c_n1, iwork, &c_n1, info);

    lwmin  = (magma_int_t)MAGMA_Z_REAL(work[0]);
    lrwmin = (magma_int_t)rwork[0];
    liwmin = (magma_int_t)iwork[0];

    if ((lwork < lwmin) && !lquery) {
        *info = -8;
    } else if ((lrwork < lrwmin) && ! lquery) {
        *info = -10;
    } else if ((liwork < liwmin) && ! lquery) {
        *info = -12;
    }

    if (*info != 0) {
	magma_xerbla("magma_zheevd", info);
	return MAGMA_ERR_ILLEGAL_VALUE;
    } else if (lquery) {
	return MAGMA_SUCCESS;
    }

    /* Quick return if possible */
    if (n == 0) {
	return MAGMA_SUCCESS;
    }

    if (n == 1) {
	w[0] = MAGMA_Z_REAL(a[0]);
	if (wantz) {
	    a[0] = MAGMA_Z_ONE;
	}
	return MAGMA_SUCCESS;
    }

    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --work;
    --rwork;
    --iwork;

    /* Get machine constants. */
    safmin = lapackf77_dlamch("Safe minimum");
    eps = lapackf77_dlamch("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = magma_dsqrt(smlnum);
    rmax = magma_dsqrt(bignum);

    /* Scale matrix to allowable range, if necessary. */
    anrm = lapackf77_zlanhe("M", uplo_, &n, &a[a_offset], &lda, &rwork[1]);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	lapackf77_zlascl(uplo_, &c__0, &c__0, &c_b18, &sigma, &n, &n, &a[a_offset], 
		&lda, info);
    }

    /* Call ZHETRD to reduce Hermitian matrix to tridiagonal form. */
    inde = 1;
    indtau = 1;
    indwrk = indtau + n;
    indrwk = inde + n;
    indwk2 = indwrk + n * n;
    llwork = lwork - indwrk + 1;
    llwrk2 = lwork - indwk2 + 1;
    llrwk = lrwork - indrwk + 1;
    /*
    lapackf77_zhetrd(uplo_, &n, &a[a_offset], &lda, &w[1], &rwork[inde], 
		     &work[indtau], &work[indwrk], &llwork, &iinfo);
    */
    magma_zhetrd(uplo_[0], n, &a[a_offset], lda, &w[1], &rwork[inde],
		 &work[indtau], &work[indwrk], llwork, &iinfo);
    
    /* For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the   
       tridiagonal matrix, then call ZUNMTR to multiply it to the Householder 
       transformations represented as Householder vectors in A. */
    if (! wantz) {
	lapackf77_dsterf(&n, &w[1], &rwork[inde], info);
    } else {
	lapackf77_zstedc("I", &n, &w[1], &rwork[inde], &work[indwrk], &n, &work[indwk2], 
		&llwrk2, &rwork[indrwk], &llrwk, &iwork[1], &liwork, info);
       
	lapackf77_zunmtr("L", uplo_, "N", &n, &n, &a[a_offset], &lda, &work[indtau], 
		&work[indwrk], &n, &work[indwk2], &llwrk2, &iinfo);
	/*
	  magma_zunmtr(MagmaLeft, uplo, MagmaNoTrans, &n, &n, &a[a_offset], &lda, &work[indtau],
	               &work[indwrk], &n, &work[indwk2], &llwrk2, &iinfo);
	*/
	lapackf77_zlacpy("A", &n, &n, &work[indwrk], &n, &a[a_offset], &lda);
    }

    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1) {
	if (*info == 0) {
	    imax = n;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	blasf77_dscal(&imax, &d__1, &w[1], &c__1);
    }

    work[1]  = MAGMA_Z_MAKE((double) lwmin, 0.);
    rwork[1] = (double) lrwmin;
    iwork[1] = liwmin;

    return MAGMA_SUCCESS;
} /* magma_zheevd */
