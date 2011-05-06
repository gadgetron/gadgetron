/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated c

*/
#include "common_magma.h"
#include <cblas.h>

#define PRECISION_c

/*
 * VERSION1 - LAPACK
 * VERSION2 - MAGMA whithout T
 * VERSION3 - MAGMA with T
 */
#define VERSION3

extern "C" magma_int_t
magma_cgeev(char jobvl, char jobvr, magma_int_t n,
	    cuFloatComplex *a, magma_int_t lda,
	    cuFloatComplex *geev_w_array,
	    cuFloatComplex *vl, magma_int_t ldvl,
	    cuFloatComplex *vr, magma_int_t ldvr,
	    cuFloatComplex *work, magma_int_t lwork,
	    float *rwork, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose   
    =======   
    CGEEV computes for an N-by-N complex nonsymmetric matrix A, the   
    eigenvalues and, optionally, the left and/or right eigenvectors.   

    The right eigenvector v(j) of A satisfies   
                     A * v(j) = lambda(j) * v(j)   
    where lambda(j) is its eigenvalue.   
    The left eigenvector u(j) of A satisfies   
                  u(j)\*\*H * A = lambda(j) * u(j)\*\*H   
    where u(j)\*\*H denotes the conjfugate transpose of u(j).   

    The computed eigenvectors are normalized to have Euclidean norm   
    equal to 1 and largest component real.   

    Arguments   
    =========   
    JOBVL   (input) CHARACTER*1   
            = 'N': left eigenvectors of A are not computed;   
            = 'V': left eigenvectors of are computed.   

    JOBVR   (input) CHARACTER*1   
            = 'N': right eigenvectors of A are not computed;   
            = 'V': right eigenvectors of A are computed.   

    N       (input) INTEGER   
            The order of the matrix A. N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.   
            On exit, A has been overwritten.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    W       (output) COMPLEX*16 array, dimension (N)   
            W contains the computed eigenvalues.   

    VL      (output) COMPLEX*16 array, dimension (LDVL,N)   
            If JOBVL = 'V', the left eigenvectors u(j) are stored one   
            after another in the columns of VL, in the same order   
            as their eigenvalues.   
            If JOBVL = 'N', VL is not referenced.   
            u(j) = VL(:,j), the j-th column of VL.   

    LDVL    (input) INTEGER   
            The leading dimension of the array VL.  LDVL >= 1; if   
            JOBVL = 'V', LDVL >= N.   

    VR      (output) COMPLEX*16 array, dimension (LDVR,N)   
            If JOBVR = 'V', the right eigenvectors v(j) are stored one   
            after another in the columns of VR, in the same order   
            as their eigenvalues.   
            If JOBVR = 'N', VR is not referenced.   
            v(j) = VR(:,j), the j-th column of VR.   

    LDVR    (input) INTEGER   
            The leading dimension of the array VR.  LDVR >= 1; if   
            JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,2*N).   
            For good performance, LWORK must generally be larger.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, the QR algorithm failed to compute all the   
                  eigenvalues, and no eigenvectors have been computed;   
                  elements and i+1:N of W contain eigenvalues which have   
                  converged.   
    =====================================================================    */


    magma_int_t c__1 = 1;
    magma_int_t c__0 = 0;
    magma_int_t c_n1 = -1;
    
    magma_int_t a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    float d__1, d__2;
    cuFloatComplex z__1, z__2;

    magma_int_t i__, k, ihi;
    float scl;
    magma_int_t ilo;
    float dum[1], eps;
    cuFloatComplex tmp;
    magma_int_t ibal;
    float anrm;
    magma_int_t ierr, itau, iwrk, nout;
    magma_int_t scalea;
    float cscale;
    magma_int_t select[1];
    static float bignum;
    static magma_int_t maxwrk;
    magma_int_t wantvl;
    static float smlnum;
    static magma_int_t irwork;
    magma_int_t lquery, wantvr;
    magma_int_t nb = 0;
    cuFloatComplex *dT = NULL;

    TimeStruct start, end;

    char side[2]   = {0, 0};
    char jobvl_[2] = {jobvl, 0};
    char jobvr_[2] = {jobvr, 0};

    *info = 0;
    lquery = lwork == -1;
    wantvl = lapackf77_lsame(jobvl_, "V");
    wantvr = lapackf77_lsame(jobvr_, "V");
    if (! wantvl && ! lapackf77_lsame(jobvl_, "N")) {
	*info = -1;
    } else if (! wantvr && ! lapackf77_lsame(jobvr_, "N")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max(1,n)) {
	*info = -5;
    } else if ( (ldvl < 1) || (wantvl && (ldvl < n))) {
	*info = -8;
    } else if ( (ldvr < 1) || (wantvr && (ldvr < n))) {
	*info = -10;
    }

    /*  Compute workspace   */
    lapackf77_cgeev(jobvl_, jobvr_, &n, a, &lda, geev_w_array, 
                    vl, &ldvl, vr, &ldvr, work, &c_n1, rwork, info);


    maxwrk = (magma_int_t)MAGMA_C_REAL(work[0]);

    if (lwork < maxwrk && ! lquery) {
        *info = -12;
    }

    if (*info != 0) {
	i__1 = -(*info);
	magma_xerbla("CGEEV ", &i__1);
	return MAGMA_ERR_ILLEGAL_VALUE;
    } else if (lquery) {
	return MAGMA_SUCCESS;
    }

    /* Quick return if possible */
    if (n == 0) {
	return MAGMA_SUCCESS;
    }
   
    // if eigenvectors are needed
#if defined(VERSION3)
    nb = magma_get_cgehrd_nb(n);
    if (CUBLAS_STATUS_SUCCESS != 
        cublasAlloc( nb*n, sizeof(cuFloatComplex), (void**)&dT)) {
	*info = -6;
	return MAGMA_ERR_CUBLASALLOC;
    }
#endif

    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    vl_dim1 = ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;

    /* Get machine constants */
    eps    = lapackf77_slamch("P");
    smlnum = lapackf77_slamch("S");
    bignum = 1. / smlnum;
    lapackf77_slabad(&smlnum, &bignum);
    smlnum = magma_ssqrt(smlnum) / eps;
    bignum = 1. / smlnum;

    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = lapackf77_clange("M", &n, &n, &a[a_offset], &lda, dum);
    scalea = 0;
    if (anrm > 0. && anrm < smlnum) {
	scalea = 1;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = 1;
	cscale = bignum;
    }
    if (scalea) {
	lapackf77_clascl("G", &c__0, &c__0, &anrm, &cscale, &n, &n, &a[a_offset], &lda, &
		ierr);
    }

    /* Balance the matrix   
       (CWorkspace: none)   
       (RWorkspace: need N) */
    ibal = 1;
    lapackf77_cgebal("B", &n, &a[a_offset], &lda, &ilo, &ihi, &rwork[ibal], &ierr);

    /* Reduce to upper Hessenberg form   
       (CWorkspace: need 2*N, prefer N+N*NB)   
       (RWorkspace: none) */
    itau = 1;
    iwrk = itau + n;
    i__1 = lwork - iwrk + 1;

    start = get_current_time();
#if defined(VERSION1)
    /*
     * Version 1 - LAPACK
     */
    lapackf77_cgehrd(&n, &ilo, &ihi, &a[a_offset], &lda,
                     &work[itau], &work[iwrk], &i__1, &ierr);
    
#elif defined(VERSION2)
    /*
     *  Version 2 - LAPACK consistent HRD
     */
    magma_cgehrd2(n, ilo, ihi, &a[a_offset], lda,
                  &work[itau], &work[iwrk], &i__1, &ierr);
    
#elif defined(VERSION3)
    /*  
     * Version 3 - LAPACK consistent MAGMA HRD + matrices T stored, 
     */
    magma_cgehrd(n, ilo, ihi, &a[a_offset], lda,
                 &work[itau], &work[iwrk], i__1, dT, &ierr);
#endif
    end = get_current_time();
    printf("    Time for cgehrd = %5.2f sec\n", GetTimerValue(start,end)/1000.);

    if (wantvl) {
      /*        Want left eigenvectors   
		Copy Householder vectors to VL */
	*(unsigned char *)side = 'L';
	lapackf77_clacpy(MagmaLowerStr, &n, &n, 
                         &a[a_offset], &lda, &vl[vl_offset], &ldvl);

	/* Generate unitary matrix in VL   
           (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)   
           (RWorkspace: none) */
	i__1 = lwork - iwrk + 1;

	start = get_current_time();
#if defined(VERSION1) || defined(VERSION2)
	/*
         * Version 1 & 2 - LAPACK
         */
        lapackf77_cunghr(&n, &ilo, &ihi, &vl[vl_offset], &ldvl, 
                         &work[itau], &work[iwrk], &i__1, &ierr);
#elif defined(VERSION3)
        /*
         * Version 3 - LAPACK consistent MAGMA HRD + matrices T stored
         */
	magma_cunghr(n, ilo, ihi, &vl[vl_offset], ldvl, &work[itau], 
		     dT, nb, &ierr);
#endif
	end = get_current_time();
	printf("    Time for cunghr = %5.2f sec\n", GetTimerValue(start,end)/1000.);

	/* Perform QR iteration, accumulating Schur vectors in VL   
           (CWorkspace: need 1, prefer HSWORK (see comments) )   
           (RWorkspace: none) */
	iwrk = itau;
	i__1 = lwork - iwrk + 1;
	lapackf77_chseqr("S", "V", &n, &ilo, &ihi, &a[a_offset], &lda, geev_w_array, &vl[
		vl_offset], &ldvl, &work[iwrk], &i__1, info);

	if (wantvr) {
	  /* Want left and right eigenvectors   
             Copy Schur vectors to VR */
	    side[0] = 'B';
	    lapackf77_clacpy("F", &n, &n, &vl[vl_offset], &ldvl, &vr[vr_offset], &ldvr);
	}

    } else if (wantvr) {
        /*  Want right eigenvectors   
            Copy Householder vectors to VR */
	side[0] = 'R';
	lapackf77_clacpy("L", &n, &n, &a[a_offset], &lda, &vr[vr_offset], &ldvr);

	/* Generate unitary matrix in VR   
           (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)   
           (RWorkspace: none) */
	i__1 = lwork - iwrk + 1;
	start = get_current_time();
#if defined(VERSION1) || defined(VERSION2)
	/*
         * Version 1 & 2 - LAPACK
         */
        lapackf77_cunghr(&n, &ilo, &ihi, &vr[vr_offset], &ldvr, 
                         &work[itau], &work[iwrk], &i__1, &ierr);
#elif defined(VERSION3)
        /*
         * Version 3 - LAPACK consistent MAGMA HRD + matrices T stored
         */
        magma_cunghr(n, ilo, ihi, &vr[vr_offset], ldvr, 
                     &work[itau], dT, nb, &ierr);
#endif
	end = get_current_time();
	printf("    Time for cunghr = %5.2f sec\n", GetTimerValue(start,end)/1000.);

	/* Perform QR iteration, accumulating Schur vectors in VR   
           (CWorkspace: need 1, prefer HSWORK (see comments) )   
           (RWorkspace: none) */
	iwrk = itau;
	i__1 = lwork - iwrk + 1;
	lapackf77_chseqr("S", "V", &n, &ilo, &ihi, &a[a_offset], &lda, geev_w_array, 
		&vr[vr_offset], &ldvr, &work[iwrk], &i__1, info);
    } else {
      /*  Compute eigenvalues only   
          (CWorkspace: need 1, prefer HSWORK (see comments) )   
          (RWorkspace: none) */
	iwrk = itau;
	i__1 = lwork - iwrk + 1;
	lapackf77_chseqr("E", "N", &n, &ilo, &ihi, &a[a_offset], &lda, geev_w_array,
		&vr[vr_offset], &ldvr, &work[iwrk], &i__1, info);
    }

    /* If INFO > 0 from ZHSEQR, then quit */
    if (*info > 0) {
	goto L50;
    }

    if (wantvl || wantvr) {
        /*  Compute left and/or right eigenvectors   
            (CWorkspace: need 2*N)   
            (RWorkspace: need 2*N) */
	irwork = ibal + n;
	lapackf77_ctrevc(side, "B", select, &n, &a[a_offset], &lda, &vl[vl_offset], &ldvl,
		&vr[vr_offset], &ldvr, &n, &nout, &work[iwrk], &rwork[irwork], 
		&ierr);
    }

    if (wantvl) {
        /*  Undo balancing of left eigenvectors   
            (CWorkspace: none)   
            (RWorkspace: need N) */
	lapackf77_cgebak("B", "L", &n, &ilo, &ihi, &rwork[ibal], &n, 
                         &vl[vl_offset], &ldvl, &ierr);

	/* Normalize left eigenvectors and make largest component real */
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scl = 1. / cblas_scnrm2(n, &vl[i__ * vl_dim1 + 1], 1);
	    cblas_csscal(n, scl, &vl[i__ * vl_dim1 + 1], 1);
	    i__2 = n;
	    for (k = 1; k <= i__2; ++k) 
	    {
                i__3 = k + i__ * vl_dim1;
                /* Computing 2nd power */
		d__1 = MAGMA_C_REAL(vl[i__3]);
		/* Computing 2nd power */
		d__2 = MAGMA_C_IMAG(vl[k + i__ * vl_dim1]);
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
            }
	    k = cblas_isamax(n, &rwork[irwork], 1);
	    MAGMA_C_CNJG(z__2, vl[k + i__ * vl_dim1]);
	    d__1 = magma_ssqrt(rwork[irwork + k - 1]);
	    MAGMA_C_DSCALE(z__1, z__2, d__1);
	    MAGMA_C_ASSIGN(tmp, z__1);
	    cblas_cscal(n, CBLAS_SADDR(tmp), &vl[i__ * vl_dim1 + 1], 1);
	    i__2 = k + i__ * vl_dim1;
	    i__3 = k + i__ * vl_dim1;
	    d__1 = MAGMA_C_REAL(vl[i__3]);
	    MAGMA_C_SET2REAL(z__1, d__1);
	    MAGMA_C_ASSIGN(vl[i__2], z__1);
	}
    }

    if (wantvr) {
      /*  Undo balancing of right eigenvectors   
          (CWorkspace: none)   
          (RWorkspace: need N) */
	lapackf77_cgebak("B", "R", &n, &ilo, &ihi, &rwork[ibal], &n, 
                         &vr[vr_offset], &ldvr, &ierr);

	/* Normalize right eigenvectors and make largest component real */
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scl = 1. / cblas_scnrm2(n, &vr[i__ * vr_dim1 + 1], 1);
	    cblas_csscal(n, scl, &vr[i__ * vr_dim1 + 1], 1);
	    i__2 = n;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = k + i__ * vr_dim1;
		/* Computing 2nd power */
		d__1 = MAGMA_C_REAL(vr[i__3]);
		/* Computing 2nd power */
		d__2 = MAGMA_C_IMAG(vr[k + i__ * vr_dim1]);
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
	    }
	    k = cblas_isamax(n, &rwork[irwork], 1);
	    MAGMA_C_CNJG(z__2, vr[k + i__ * vr_dim1]);
	    d__1 = magma_ssqrt(rwork[irwork + k - 1]);
	    MAGMA_C_DSCALE(z__1, z__2, d__1);
	    MAGMA_C_ASSIGN(tmp, z__1);
	    cblas_cscal(n, CBLAS_SADDR(tmp), &vr[i__ * vr_dim1 + 1], 1);
	    i__2 = k + i__ * vr_dim1;
	    i__3 = k + i__ * vr_dim1;
	    d__1 = MAGMA_C_REAL(vr[i__3]);
	    MAGMA_C_SET2REAL(z__1, d__1);
	    MAGMA_C_ASSIGN(vr[i__2], z__1);
	}
    }

    /*  Undo scaling if necessary */
L50:
    if (scalea) {
	i__1 = n - *info;
	/* Computing MAX */
	i__3 = n - *info;
	i__2 = max(i__3,1);
	lapackf77_clascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, 
		geev_w_array + *info, &i__2, &ierr);
	if (*info > 0) {
	    i__1 = ilo - 1;
	    lapackf77_clascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, 
		    geev_w_array, &n, &ierr);
	}
    }

#if defined(VERSION3)
    cublasFree( dT );
#endif
    work[1] = MAGMA_C_MAKE((float) maxwrk, 0.);
    return MAGMA_SUCCESS;
} /* magma_cgeev */

