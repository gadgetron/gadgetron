/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions mixed zc -> ds

*/
#include "common_magma.h"

#define BWDMAX 1.0
#define ITERMAX 30

extern "C" magma_int_t
magma_zcgeqrsv_gpu(magma_int_t M, magma_int_t N, magma_int_t NRHS, 
                   cuDoubleComplex *dA,  magma_int_t ldda, 
                   cuDoubleComplex *dB,  magma_int_t lddb, 
                   cuDoubleComplex *dX,  magma_int_t lddx, 
		   magma_int_t *iter, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    ZCGEQRSV solves the least squares problem 
       min || A*X - B ||,
    where A is an M-by-N matrix and X and B are M-by-NRHS matrices.

    ZCGEQRSV first attempts to factorize the matrix in SINGLE PRECISION
    and use this factorization within an iterative refinement procedure
    to produce a solution with DOUBLE PRECISION norm-wise backward error
    quality (see below). If the approach fails the method switches to a
    DOUBLE PRECISION factorization and solve.

    The iterative refinement is not going to be a winning strategy if
    the ratio SINGLE PRECISION performance over DOUBLE PRECISION
    performance is too small. A reasonable strategy should take the
    number of right-hand sides and the size of the matrix into account.
    This might be done with a call to ILAENV in the future. Up to now, we
    always try iterative refinement.
    The iterative refinement process is stopped if
        ITER > ITERMAX
    or for all the RHS we have:
        RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
    where
        o ITER is the number of the current iteration in the iterative
          refinement process
        o RNRM is the infinity-norm of the residual
        o XNRM is the infinity-norm of the solution
        o ANRM is the infinity-operator-norm of the matrix A
        o EPS is the machine epsilon returned by DLAMCH('Epsilon')
    The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 respectively.

    Arguments
    =========

    M       (input) INTEGER   
            The number of rows of the matrix A. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A. M >= N >= 0.

    NRHS    (input) INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    A       (input or input/output) DOUBLE PRECISION array, dimension (ldda,N)
            On entry, the M-by-N coefficient matrix A.
            On exit, if iterative refinement has been successfully used
            (info.EQ.0 and ITER.GE.0, see description below), A is
            unchanged. If double precision factorization has been used
            (info.EQ.0 and ITER.LT.0, see description below), then the
            array A contains the QR factorization of A as returned by
            function DGEQRF_GPU.

    ldda     (input) INTEGER
            The leading dimension of the array A.  ldda >= max(1,M).

    B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
            The M-by-NRHS right hand side matrix B.

    LDB     (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,M).

    X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)
            If info = 0, the N-by-NRHS solution matrix X.

    LDX     (input) INTEGER
            The leading dimension of the array X.  LDX >= max(1,N).

    WORK    (workspace) DOUBLE PRECISION array, dimension (N*NRHS)
            This array is used to hold the residual vectors.

    SWORK   (workspace) REAL array, dimension (M*(N+NRHS))
            This array is used to store the single precision matrix and the
            right-hand sides or solutions in single precision.

    ITER    (output) INTEGER
            < 0: iterative refinement has failed, double precision
                 factorization has been performed
                 -1 : the routine fell back to full precision for
                      implementation- or machine-specific reasons
                 -2 : narrowing the precision induced an overflow,
                      the routine fell back to full precision
                 -3 : failure of SGETRF
                 -31: stop the iterative refinement after the 30th
                      iterations
            > 0: iterative refinement has been successfully used.
                 Returns the number of iterations
 
    info    (output) INTEGER
            = 0:  successful exit
            < 0:  if info = -i, the i-th argument had an illegal value

    TAU     (output) REAL array, dimension (N)
            On exit, TAU(i) contains the scalar factor of the elementary
            reflector H(i), as returned by magma_cgeqrf_gpu.

    LWORK   (input) INTEGER   
            The dimension of the array H_WORK.  LWORK >= (M+N+NB)*NB,   
            where NB can be obtained through magma_get_sgeqrf_nb(M).

    H_WORK  (workspace/output) REAL array, dimension (MAX(1,LWORK))   
            Higher performance is achieved if H_WORK is in pinned memory, e.g.
            allocated using cudaMallocHost.

    D_WORK  (workspace/output)  REAL array on the GPU, dimension 2*N*NB,
            where NB can be obtained through magma_get_sgeqrf_nb(M).
            It starts with NB*NB blocks that store the triangular T 
            matrices, followed by the NB*NB blocks of the diagonal 
            inverses for the R matrix.

    TAU_D   (output) DOUBLE REAL array, dimension (N)
            On exit, if the matrix had to be factored in double precision,
            TAU(i) contains the scalar factor of the elementary
            reflector H(i), as returned by magma_zgeqrf_gpu.

    LWORK_D (input) INTEGER   
            The dimension of the array H_WORK_D. LWORK_D >= (M+N+NB)*NB,   
            where NB can be obtained through magma_get_dgeqrf_nb(M).

    H_WORK_D (workspace/output) DOUBLE REAL array, dimension (MAX(1,LWORK_D))
            This memory is unattached if the iterative refinement worked, 
            otherwise it is used as workspace to factor the matrix in
            double precision. Higher performance is achieved if H_WORK_D is 
            in pinned memory, e.g. allocated using cudaMallocHost. 

    D_WORK_D (workspace/output) DOUBLE REAL array on the GPU, dimension 2*N*NB,
            where NB can be obtained through magma_get_dgeqrf_nb(M).
            This memory is unattached if the iterative refinement worked, 
            otherwise it is used as workspace to factor the matrix in
            double precision. It starts with NB*NB blocks that store the 
            triangular T matrices, followed by the NB*NB blocks of the 
            diagonal inverses for the R matrix.

    =====================================================================    */

    cuDoubleComplex mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex zone  = MAGMA_Z_ONE;
    magma_int_t     ione  = 1;
    cuDoubleComplex *dworkd, *hworkd;
    cuFloatComplex  *dworks, *hworks;
    cuDoubleComplex *dR, *tau, *dT;
    cuFloatComplex  *dSA, *dSX, *dST, *stau;
    cuDoubleComplex Xnrmv, Rnrmv; 
    double          Anrm, Xnrm, Rnrm, cte, eps; 
    magma_int_t     i, j, iiter, nb, lhwork, minmn, size, ret;
    
    /*
      Check The Parameters. 
    */
    *iter = 0 ;
    *info = 0 ;
    if ( N < 0 )
        *info = -1;
    else if(NRHS<0)
        *info = -3;
    else if( ldda < max(1,N))
        *info = -5;
    else if( lddb < max(1,N))
        *info = -7;
    else if( lddx < max(1,N))
        *info = -9;

    if( *info != 0 ){
        magma_xerbla("magma_zcgeqrsv_gpu", info);
	return MAGMA_ERR_ILLEGAL_VALUE;
    }

    if( N == 0 || NRHS == 0 )
        return MAGMA_SUCCESS;

    nb   = magma_get_cgeqrf_nb(M);
    minmn= min(M, N);

    /*
     * Allocate temporary buffers
     */
    /* dworks(dSA + dSX + dST) */
    size = ldda*N +  N*NRHS 
	+  ( 2*minmn + ((N+31)/32)*32 )*nb;
    if( CUBLAS_STATUS_SUCCESS != cublasAlloc(size, sizeof(cuFloatComplex), (void**)&dworks) ) {
	fprintf(stderr, "Allocation of dworks failed (%d)\n", size);
        magma_xerbla("magma_zcgeqrsv_gpu", info);
	return MAGMA_ERR_CUBLASALLOC;
    }
    dSA = dworks;
    dSX = dSA + ldda*N;
    dST = dSX + N*NRHS;
    
    /* dworkd(dR) = N*NRHS */
    size = N*NRHS;
    if( CUBLAS_STATUS_SUCCESS != cublasAlloc(size, sizeof(cuDoubleComplex), (void**)&dworkd) ) {
	cublasFree(dworks);
	fprintf(stderr, "Allocation of dworkd failed\n");
        magma_xerbla("magma_zcgeqrsv_gpu", info);
	return MAGMA_ERR_CUBLASALLOC;
    }
    dR = dworkd;

    /* hworks(stau + workspace for cgeqrs) = min(M,N) + lhworks */
    lhwork = nb*max((M-N+nb+2*(NRHS)), 1);
    lhwork = max(lhwork, N*nb); /* We hope that magma nb is bigger than lapack nb to have enough memory in workspace */
    size = minmn + lhwork;
    hworks = (cuFloatComplex*) malloc( size * sizeof(cuFloatComplex) );
    if( hworks == NULL ) {
	cublasFree(dworks);
	cublasFree(dworkd);
	fprintf(stderr, "Allocation of hworks failed\n");
        magma_xerbla("magma_zcgeqrsv_gpu", info);
	return MAGMA_ERR_ALLOCATION;
    }
    stau = hworks;
    hworks += minmn;

    eps  = lapackf77_dlamch("Epsilon");
    Anrm = magmablas_zlange('I', M, N, dA, ldda, (double*)dworkd );
    cte  = Anrm * eps *  pow((double)N, 0.5) * BWDMAX ;

    /*
     * Convert to single precision
     */
    magmablas_zlag2c(N, NRHS, dB, lddb, dSX, N, info );
    if( *info != 0 ) {
	*iter = -2; goto L40;
    }

    magmablas_zlag2c(N, N, dA, ldda, dSA, ldda, info );
    if(*info !=0){
        *iter = -2; goto L40;
    }

    // In an ideal version these variables should come from user.
    magma_cgeqrf_gpu(M, N, dSA, ldda, stau, dST, info);
    if( *info != 0 ) {
        *iter = -3; goto L40;
    }

    magma_cgeqrs_gpu(M, N, NRHS, dSA, ldda, stau, dST, dSX, N, hworks, lhwork, info);

    // dX = dSX
    magmablas_clag2z(N, NRHS, dSX, N, dX, lddx, info);

    // dR = dB
    magmablas_zlacpy(MagmaUpperLower, N, NRHS, dB, lddb, dR, N);

    // dR = dB - dA * dX
    if( NRHS == 1 )
        cublasZgemv( MagmaNoTrans, N, N, 
		     mzone, dA, ldda, 
		            dX, 1, 
		     zone,  dR, 1);
    else
        cublasZgemm( MagmaNoTrans, MagmaNoTrans, N, NRHS, N, 
                     mzone, dA, ldda, 
                            dX, lddx, 
                     zone,  dR, N );

    for(i=0; i<NRHS; i++){
        j = cublasIzamax( N, dX+i*N, 1);
        cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dX+i*N+j-1, 1, &Xnrmv, 1);
        Xnrm = lapackf77_zlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
      
        j = cublasIzamax ( N, dR+i*N, 1 );
        cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dR+i*N+j-1, 1, &Rnrmv, 1);
        Rnrm = lapackf77_zlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
      
        if( Rnrm >  Xnrm *cte ) goto L10;
    }

    *iter = 0;

    /* Free workspaces */
    cublasFree(dworks);
    cublasFree(dworkd);
    free(stau);
    return MAGMA_SUCCESS;

  L10:
    for(iiter=1; iiter<ITERMAX; ) {
        *info = 0 ;
        /*  Convert R from double precision to single precision
            and store the result in SX.
            Solve the system SA*SX = SR.
            -- These two Tasks are merged here. */
        // make SWORK = WORK ... residuals... 
        magmablas_zlag2c( N, NRHS, dR, N, dSX, N, info );
        magma_cgeqrs_gpu( M, N, NRHS, dSA, ldda, stau, dST, dSX, N, hworks, lhwork, info);

        if( *info != 0 ){
            *iter = -3; goto L40;
        }

        for(i=0; i<NRHS; i++) {
            magmablas_zcaxpycp( dSX+i*N, dX+i*lddx, N, dB+i*lddb, dR+i*N );
        }

        /* unnecessary may be */
        magmablas_zlacpy(MagmaUpperLower, N, NRHS, dB, lddb, dR, N);
        if( NRHS == 1 )
            cublasZgemv( MagmaNoTrans, N, N, 
			 mzone, dA, ldda,
			        dX, 1,
			 zone,  dR, 1);
        else
            cublasZgemm( MagmaNoTrans, MagmaNoTrans, N, NRHS, N, 
                         mzone, dA, ldda,
			        dX, lddx,
			 zone,  dR, N);

        /*  Check whether the NRHS normwise backward errors satisfy the
            stopping criterion. If yes, set ITER=IITER>0 and return.     */
        for(i=0;i<NRHS;i++)
        {
	    j = cublasIzamax( N, dX+i*N, 1);
	    cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dX+i*N+j-1, 1, &Xnrmv, 1);
	    Xnrm = lapackf77_zlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
	    
	    j = cublasIzamax ( N, dR+i*N, 1 );
	    cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dR+i*N+j-1, 1, &Rnrmv, 1);
	    Rnrm = lapackf77_zlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
	    
	    if( Rnrm >  Xnrm *cte ) goto L20;
        }

        /*  If we are here, the NRHS normwise backward errors satisfy
            the stopping criterion, we are good to exit.                    */
        *iter = iiter ;

	/* Free workspaces */
	cublasFree(dworks);
	cublasFree(dworkd);
	free(stau);
        return MAGMA_SUCCESS;
      L20:
        iiter++;
    }

    /* If we are at this place of the code, this is because we have
       performed ITER=ITERMAX iterations and never satisified the
       stopping criterion, set up the ITER flag accordingly and follow
       up on double precision routine.                                    */
    *iter = -ITERMAX - 1 ;

  L40:
    cublasFree(dworks);

    /*
     * Allocate temporary buffers
     */
    /* dworkd(dT + tau) = min_mn + min_mn*nb*3 */
    nb   = magma_get_zgeqrf_nb(M);
    size = minmn * (3 * nb + 1);
    if ( size > (N*NRHS) ) {
	cublasFree(dworkd);
	if( CUBLAS_STATUS_SUCCESS != cublasAlloc(size, sizeof(cuDoubleComplex), (void**)&dworkd) ) {
	    fprintf(stderr, "Allocation of dworkd2 failed\n");
	    magma_xerbla("magma_zcgeqrsv_gpu", info);
	    return MAGMA_ERR_CUBLASALLOC;
	}
    }
    tau = dworkd;
    dT  = tau + minmn;

    /* hworks(stau + workspace for cgeqrs) = min(M,N) + lhworks */
    if ( (2*lhwork) > (minmn+lhwork) ) {
	free(stau);
	hworks = (cuFloatComplex*) malloc( lhwork * sizeof(cuDoubleComplex) );
	if( hworks == NULL ) {
	    cublasFree(dworkd);
	    fprintf(stderr, "Allocation of hworkd2 failed\n");
	    magma_xerbla("magma_zcgeqrsv_gpu", info);
	    return MAGMA_ERR_ALLOCATION;
	}
    }
    hworkd = (cuDoubleComplex*) hworks;

    /* Single-precision iterative refinement failed to converge to a
       satisfactory solution, so we resort to double precision.           */
    ret = magma_zgeqrf_gpu(M, N, dA, ldda, tau, dT, info);
    if( (ret != MAGMA_SUCCESS) || (*info != 0) ){
	cublasFree(dworkd);
	free(hworkd);
	return ret;
    }
    magmablas_zlacpy(MagmaUpperLower, N, NRHS, dB, lddb, dX, lddx);
    ret = magma_zgeqrs_gpu(M, N, NRHS, dA, ldda, tau, dT, dX, lddx, hworkd, lhwork, info);
    cublasFree(dworkd);
    free(hworkd);
    return ret;
}

