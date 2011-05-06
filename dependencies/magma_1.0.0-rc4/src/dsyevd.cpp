/*    
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal d -> s

*/
#include "common_magma.h"

extern "C" magma_int_t 
magma_dsyevd(char jobz, char uplo, 
	     magma_int_t n, 
	     double *a, magma_int_t lda, 
	     double *w, 
	     double *work, magma_int_t lwork,
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
    DSYEVD computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric matrix A.  If eigenvectors are desired, it uses a   
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

    A       (input/output) DOUBLE_PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
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

    WORK    (workspace/output) DOUBLE_PRECISION array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.   
            If N <= 1,                LWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LWORK must be at least 2*N + 1.   
            If JOBZ  = 'V' and N > 1, LWORK must be at least 1 + 6*N + 2*N**2.   

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal sizes of the WORK and IWORK
            arrays, returns these values as the first entries of the WORK
            and IWORK arrays, and no error message related to LWORK or
            LIWORK is issued by XERBLA.

    IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))   
            On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.   

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If N <= 1,                LIWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.   

            If LIWORK = -1, then a workspace query is assumed; the
            routine only calculates the optimal sizes of the WORK and
            IWORK arrays, returns these values as the first entries of
            the WORK and IWORK arrays, and no error message related to
            LWORK or LIWORK is issued by XERBLA.

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
    static double rmin, rmax;
    static magma_int_t lopt;
    static double sigma;
    static magma_int_t iinfo, lwmin, liopt;
    static magma_int_t lower;
    static magma_int_t wantz;
    static magma_int_t indwk2, llwrk2;
    static magma_int_t iscale;
    static double safmin;
    static double bignum;
    static magma_int_t indtau;
    static magma_int_t indwrk, liwmin;
    static magma_int_t llwork;
    static double smlnum;
    static magma_int_t lquery;

    /* Function Body */
    wantz = lapackf77_lsame(jobz_, MagmaVectorsStr);
    lower = lapackf77_lsame(uplo_, MagmaLowerStr);
    lquery = lwork == -1 || liwork == -1;

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

    lapackf77_dsyevd(jobz_, uplo_, &n, 
                     a, &lda, w, work, &c_n1, 
                     iwork, &c_n1, info);

    lwmin  = (magma_int_t)work[0];
    liwmin = (magma_int_t)iwork[0];

    if ((lwork < lwmin) && !lquery) {
        *info = -8;
    } else if ((liwork < liwmin) && ! lquery) {
        *info = -10;
    }

    if (*info != 0) {
	magma_xerbla("magma_dsyevd", info);
	return MAGMA_ERR_ILLEGAL_VALUE;
    } else if (lquery) {
	return MAGMA_SUCCESS;
    }

    /* Quick return if possible */
    if (n == 0) {
	return MAGMA_SUCCESS;
    }

    if (n == 1) {
        w[0] = a[0];
	if (wantz) {
            a[0] = 1.;
	}
	return MAGMA_SUCCESS;
    }

    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --work;
    --iwork;

    /* Get machine constants. */
    safmin = lapackf77_dlamch("Safe minimum");
    eps = lapackf77_dlamch("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = magma_dsqrt(smlnum);
    rmax = magma_dsqrt(bignum);

    /* Scale matrix to allowable range, if necessary. */
    anrm = lapackf77_dlansy("M", uplo_, &n, &a[a_offset], &lda, &work[1]);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	lapackf77_dlascl(uplo_, &c__0, &c__0, &c_b18, &sigma, &n, &n, &a[a_offset], 
		&lda, info);
    }

    /* Call DSYTRD to reduce symmetric matrix to tridiagonal form. */
    inde = 1;
    indtau = inde + n;
    indwrk = indtau + n;
    llwork = lwork - indwrk + 1;
    indwk2 = indwrk + n * n;
    llwrk2 = lwork - indwk2 + 1;
  
    /*
    lapackf77_dsytrd(uplo_, &n, &a[a_offset], &lda, &w[1], &work[inde], 
		     &work[indtau], &work[indwrk], &llwork, &iinfo);
    */
    magma_dsytrd(uplo_[0], n, &a[a_offset], lda, &w[1], &work[inde],
		 &work[indtau], &work[indwrk], llwork, &iinfo);
    
    /* For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the   
       tridiagonal matrix, then call DORMTR to multiply it to the Householder 
       transformations represented as Householder vectors in A. */
    if (! wantz) {
	lapackf77_dsterf(&n, &w[1], &work[inde], info);
    } else {
	lapackf77_dstedc("I", &n, &w[1], &work[inde], &work[indwrk], &n, &work[indwk2], 
		&llwrk2, &iwork[1], &liwork, info);
       
	lapackf77_dormtr("L", uplo_, "N", &n, &n, &a[a_offset], &lda, &work[indtau], 
		&work[indwrk], &n, &work[indwk2], &llwrk2, &iinfo);
	/*
	  magma_dormtr("L", uplo_, "N", n, n, &a[a_offset], &lda, &work[indtau],
	               &work[indwrk], n, &work[indwk2], &llwrk2, &iinfo);
	*/
	lapackf77_dlacpy("A", &n, &n, &work[indwrk], &n, &a[a_offset], &lda);
    }

    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1) {
	d__1 = 1. / sigma;
	blasf77_dscal(&n, &d__1, &w[1], &c__1);
    }

    MAGMA_D_SET2REAL(work[1], (double) lopt);
    iwork[1] = liopt;

    return MAGMA_SUCCESS;
} /* magma_dsyevd */
