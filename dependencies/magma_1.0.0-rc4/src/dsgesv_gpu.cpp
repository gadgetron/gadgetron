/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated ds

*/
#include "common_magma.h"

#define BWDMAX 1.0
#define ITERMAX 30

extern "C" magma_int_t
magma_dsgesv_gpu(char trans, magma_int_t N, magma_int_t NRHS, 
                 double *dA, magma_int_t ldda, 
                 magma_int_t *IPIV,  magma_int_t *dIPIV,
                 double *dB, magma_int_t lddb, 
                 double *dX, magma_int_t lddx, 
                 double *dworkd, float *dworks,
		 magma_int_t *iter, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    DSGESV computes the solution to a complex system of linear equations
       A * X = B or A' * X = B
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

    DSGESV first attempts to factorize the matrix in SINGLE PRECISION
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
    TRANS   (input) CHARACTER*1
            Specifies the form of the system of equations:
            = 'N':  A * X = B  (No transpose)
            = 'T':  A'* X = B  (Transpose)
            = 'C':  A'* X = B  (ugate transpose = Transpose)

    N       (input) INTEGER
            The number of linear equations, i.e., the order of the
            matrix A.  N >= 0.

    NRHS    (input) INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    dA      (input or input/output) DOUBLE PRECISION array, dimension (ldda,N)
            On entry, the N-by-N coefficient matrix A.
            On exit, if iterative refinement has been successfully used
            (info.EQ.0 and ITER.GE.0, see description below), A is
            unchanged. If double precision factorization has been used
            (info.EQ.0 and ITER.LT.0, see description below), then the
            array A contains the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    ldda    (input) INTEGER
            The leading dimension of the array A.  ldda >= max(1,N).

    IPIV    (output) INTEGER array, dimension (N)
            The pivot indices that define the permutation matrix P;
            row i of the matrix was interchanged with row IPIV(i).
            Corresponds either to the single precision factorization
            (if info.EQ.0 and ITER.GE.0) or the double precision
            factorization (if info.EQ.0 and ITER.LT.0).

    dIPIV   (output) INTEGER array on the GPU, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was moved to row IPIV(i).

    dB      (input) DOUBLE PRECISION array, dimension (lddb,NRHS)
            The N-by-NRHS right hand side matrix B.

    lddb    (input) INTEGER
            The leading dimension of the array B.  lddb >= max(1,N).

    dX      (output) DOUBLE PRECISION array, dimension (lddx,NRHS)
            If info = 0, the N-by-NRHS solution matrix X.

    lddx    (input) INTEGER
            The leading dimension of the array X.  lddx >= max(1,N).

    dworkd  (workspace) DOUBLE PRECISION array, dimension (N*NRHS)
            This array is used to hold the residual vectors.

    dworks  (workspace) SINGLE PRECISION array, dimension (N*(N+NRHS))
            This array is used to store the single precision matrix and the
            right-hand sides or solutions in single precision.

    iter    (output) INTEGER
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
 
    info   (output) INTEGER
            = 0:  successful exit
            < 0:  if info = -i, the i-th argument had an illegal value
            > 0:  if info = i, U(i,i) computed in DOUBLE PRECISION is
                  exactly zero.  The factorization has been completed,
                  but the factor U is exactly singular, so the solution
                  could not be computed.
    =====================================================================    */

    double mzone = MAGMA_D_NEG_ONE;
    double zone  = MAGMA_D_ONE;
    magma_int_t     ione  = 1;
    double          cte, eps;
    double Xnrmv, Rnrmv;
    float  *dSA, *dSX;
    double          Anrm, Xnrm, Rnrm;
    magma_int_t     i, j, iiter, ret;
    
    /*
      Check The Parameters. 
    */
    *info = *iter = 0 ;
    if ( N <0)
	*info = -1;
    else if(NRHS<0)
	*info =-2;
    else if(ldda < max(1,N))
	*info =-4;
    else if( lddb < max(1,N))
	*info =-8;
    else if( lddx < max(1,N))
	*info =-10;
    
    if(*info!=0){
	magma_xerbla("magma_dsgesv_gpu",info);
	return MAGMA_ERR_ILLEGAL_VALUE;
    }
    
    if( N == 0 || NRHS == 0 )
	return MAGMA_SUCCESS;
    
    eps  = lapackf77_dlamch("Epsilon");
    Anrm = magmablas_dlange('I', N, N, dA, ldda, (double*)dworkd );
    cte  = Anrm * eps * pow((double)N,0.5) * BWDMAX;
    
    dSX = dworks;
    dSA = dworks+N*NRHS;
    
    if(*info !=0){
	*iter = -2 ;
	printf("magmablas_dlag2s\n");
	goto L40;
    }
    magmablas_dlag2s(N, N, dA, ldda, dSA, N, info ); // Merge with DLANGE /
    if(*info !=0){
	*iter = -2 ;
	printf("magmablas_dlag2s\n");
	goto L40;
    }
    
    magma_sgetrf_gpu(N, N, dSA, N, IPIV, info);
    
    // Generate parallel pivots
    {
	int *newipiv  = (int*)malloc(N * sizeof(int));
	swp2pswp(trans, N, IPIV, newipiv);
	cudaMemcpy(dIPIV, newipiv, N*sizeof(int), cudaMemcpyHostToDevice);
	free(newipiv);
    }
    
    if(info[0] !=0){
	*iter = -3 ;
	goto L40;
    }
    magma_dsgetrs_gpu(trans, N, NRHS, dSA, N, dIPIV, dB, lddb, dX, lddx, dSX, info);
    
    magmablas_dlacpy(MagmaUpperLower, N, NRHS, dB, lddb, dworkd, N);
     /* TODO: update optimisation from gemv_MLU into classic gemv */
    if ( NRHS == 1 ) 
	cublasDgemv( trans, N, N, mzone, dA, ldda, dX, 1, zone, dworkd, 1);
    else
	cublasDgemm( trans, MagmaNoTrans, N, NRHS, N, mzone, dA, ldda, dX, lddx, zone, dworkd, N);
    
    for(i=0;i<NRHS;i++)
    {
	j = cublasIdamax( N, dX+i*N, 1) ;
	cublasGetMatrix( 1, 1, sizeof(double), dX+i*N+j-1, 1, &Xnrmv, 1);
	Xnrm = lapackf77_dlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
	
	j = cublasIdamax ( N, dworkd+i*N, 1 );
	cublasGetMatrix( 1, 1, sizeof(double), dworkd+i*N+j-1, 1, &Rnrmv, 1 );
	Rnrm = lapackf77_dlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
	
	if( Rnrm >  (Xnrm*cte) ){
	    goto L10;
	}
    }
    
    *iter = 0;
    return MAGMA_SUCCESS;
  L10:
    ;
    
    for(iiter=1;iiter<ITERMAX;)
    {
	*info = 0 ;
	/*
	  Convert R (in dworkd) from double precision to single precision
	  and store the result in SX.
	  Solve the system SA * X = dworkd and store the result in dworkd.
	  -- These two Tasks are merged here. 
	*/
	magma_dsgetrs_gpu(trans, N, NRHS, dSA, N, dIPIV, dworkd, lddb, dworkd, N, dSX, info);
	if(info[0] != 0){
	    *iter = -3 ;
	    goto L40;
	}
	/* Add the correction (dworkd) to dX and make dworkd = dB.
           This saves going through dworkd a second time (if done with one more kernel). */
	for(i=0;i<NRHS;i++){
	    magmablas_daxpycp(dworkd+i*N, dX+i*N, N, dB+i*N);
	}
	
	//magmablas_dlacpy(MagmaUpperLower, N, NRHS, dB, lddb, dworkd, N);
	if( NRHS == 1 )
	    /* TODO: update optimisation from gemv_MLU into classic gemv */
	    cublasDgemv( trans, N, N, mzone, dA, ldda, dX, 1, zone, dworkd, 1);
	else
	    cublasDgemm( trans, MagmaNoTrans, N, NRHS, N, mzone, dA, ldda, dX, lddx, zone, dworkd, N);
	
	/*
	  Check whether the NRHS normwise backward errors satisfy the
	  stopping criterion. If yes, set ITER=IITER>0 and return.
	*/
	for(i=0;i<NRHS;i++)
	{
	    j = cublasIdamax( N, dX+i*N, 1) ;
	    cublasGetMatrix( 1, 1, sizeof(double), dX+i*N+j-1, 1, &Xnrmv, 1);
	    Xnrm = lapackf77_dlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
	    
	    j = cublasIdamax ( N, dworkd+i*N, 1 );
	    cublasGetMatrix( 1, 1, sizeof(double), dworkd+i*N+j-1, 1, &Rnrmv, 1 );
	    Rnrm = lapackf77_dlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
	    
	    if( Rnrm >  Xnrm *cte ){
		goto L20;
	    }
	}
	/*
	  If we are here, the NRHS normwise backward errors satisfy the
	  stopping criterion, we are good to exit.
	*/
	
	*iter = iiter ;
	return MAGMA_SUCCESS;
      L20:
	iiter++ ;
    }
    /*
      If we are at this place of the code, this is because we have
      performed ITER=ITERMAX iterations and never satisified the
      stopping criterion, set up the ITER flag accordingly and follow up
      on double precision routine.
    */
    *iter = -ITERMAX - 1 ;
    
  L40:
    /*
      Single-precision iterative refinement failed to converge to a
      satisfactory solution, so we resort to double precision.  
    */
    if( *info != 0 ){
	return MAGMA_SUCCESS;
    }
    
    ret = magma_dgetrf_gpu( N, N, dA, ldda, IPIV, info );
    if( (ret != MAGMA_SUCCESS) || (*info != 0) ){
	return ret;
    }
    magmablas_dlacpy( MagmaUpperLower, N, NRHS, dB, lddb, dX, N );
    ret = magma_dgetrs_gpu( trans, N, NRHS, dA, ldda, IPIV, dX, N, info );
    return ret;
}

