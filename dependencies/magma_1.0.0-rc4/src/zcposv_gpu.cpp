/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions mixed zc -> ds

*/
#include "common_magma.h"

#define ITERMAX 30
#define BWDMAX 1.0

// === Define what BLAS to use ============================================
#define cublasZhemv magmablas_zhemv
// === End defining what BLAS to use ======================================

extern "C" magma_int_t
magma_zcposv_gpu(char uplo, magma_int_t n, magma_int_t nrhs, 
                 cuDoubleComplex *dA, magma_int_t ldda, 
                 cuDoubleComplex *dB, magma_int_t lddb, 
                 cuDoubleComplex *dX, magma_int_t lddx, 
                 cuDoubleComplex *dworkd, cuFloatComplex *dworks,
                 magma_int_t *iter, magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    ZCPOSV computes the solution to a complex system of linear equations
       A * X = B,
    where A is an N-by-N symmetric positive definite matrix and X and B
    are N-by-NRHS matrices.

    ZCPOSV first attempts to factorize the matrix in complex SINGLE PRECISION
    and use this factorization within an iterative refinement procedure
    to produce a solution with complex DOUBLE PRECISION norm-wise backward error
    quality (see below). If the approach fails the method switches to a
    complex DOUBLE PRECISION factorization and solve.

    The iterative refinement is not going to be a winning strategy if
    the ratio complex SINGLE PRECISION performance over DOUBLE PRECISION
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

    UPLO    (input) CHARACTER
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.

    N       (input) INTEGER
            The number of linear equations, i.e., the order of the
            matrix A.  N >= 0.

    NRHS    (input) INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    A       (input or input/output) COMPLEX_16 array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = 'L', the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit, if iterative refinement has been successfully used
            (INFO.EQ.0 and ITER.GE.0, see description below), then A is
            unchanged, if double factorization has been used
            (INFO.EQ.0 and ITER.LT.0, see description below), then the
            array A contains the factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    B       (input) COMPLEX_16 array, dimension (LDB,NRHS)
            The N-by-NRHS right hand side matrix B.

    LDB     (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    X       (output) COMPLEX_16 array, dimension (LDX,NRHS)
            If INFO = 0, the N-by-NRHS solution matrix X.

    LDX     (input) INTEGER
            The leading dimension of the array X.  LDX >= max(1,N).

    dworkd  (workspace) COMPLEX_16 array, dimension (N*NRHS)
            This array is used to hold the residual vectors.

    dworks  (workspace) COMPLEX array, dimension (N*(N+NRHS))
            This array is used to use the complex single precision matrix 
            and the right-hand sides or solutions in single precision.

    ITER    (output) INTEGER
            < 0: iterative refinement has failed, double precision
                 factorization has been performed
                 -1 : the routine fell back to full precision for
                      implementation- or machine-specific reasons
                 -2 : narrowing the precision induced an overflow,
                      the routine fell back to full precision
                 -3 : failure of SPOTRF
                 -31: stop the iterative refinement after the 30th
                      iterations
            > 0: iterative refinement has been successfully used.
                 Returns the number of iterations

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  if INFO = i, the leading minor of order i of (DOUBLE
                  PRECISION) A is not positive definite, so the
                  factorization could not be completed, and the solution
                  has not been computed.

    =====================================================================    */


    cuDoubleComplex mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex zone  = MAGMA_Z_ONE;
    magma_int_t     ione  = 1;
    cuFloatComplex *dSA, *dSX;
    cuDoubleComplex Xnrmv, Rnrmv; 
    double          Xnrm, Rnrm, Anrm, cte, eps;
    magma_int_t     i, j, iiter, ret;

    *iter = 0 ;
    *info = 0 ; 

    if ( n <0)
        *info = -1;
    else if( nrhs<0 )
        *info =-2;
    else if( ldda < max(1,n) )
        *info =-4;
    else if( lddb < max(1,n) )
        *info =-7;
    else if( lddx < max(1,n) )
        *info =-9;
   
    if( *info != 0 ){
        magma_xerbla("magma_zcposv", info);	
	return MAGMA_ERR_ILLEGAL_VALUE;
    }

    if( n == 0 || nrhs == 0 ) 
        return MAGMA_SUCCESS;

    eps = lapackf77_dlamch("Epsilon");

    //ANRM = magmablas_zlanhe(  'I',  uplo , N ,A, LDA, (double*)dworkd);
    //cte  = ANRM * EPS *  pow((double)N,0.5) * BWDMAX ;  

    dSX = dworks;
    dSA = dworks + n*nrhs;
 
    magmablas_zlag2c(n, nrhs, dB, lddb, dSX, n, info );
    if( *info !=0 ){
        *iter = -2;
        goto L40;
    } 
  
    magmablas_zlat2c(uplo, n, dA, ldda, dSA, n, info ); 
    if( *info !=0 ){
        *iter = -2 ;
        goto L40;
    }
 
    Anrm = magmablas_clanhe( 'I', uplo, n, dSA, n, (float *)dworkd);
    cte  = Anrm * eps * pow((double)n,0.5) * BWDMAX ;

    magma_cpotrf_gpu(uplo, n, dSA, ldda, info);
    if( *info !=0 ){
        *iter = -3 ;
        goto L40;
    }
    magma_cpotrs_gpu(uplo, n, nrhs, dSA, ldda, dSX, lddb, info);
    magmablas_clag2z(n, nrhs, dSX, n, dX, lddx, info);
  
    magmablas_zlacpy( MagmaUpperLower, n, nrhs, dB, lddb, dworkd, n);
    
    if( nrhs == 1 )
        cublasZhemv(uplo, n, mzone, dA, ldda, dX, 1, zone, dworkd, 1);
    else
        cublasZhemm(MagmaLeft, uplo, n, nrhs, mzone, dA, ldda, dX, lddx, zone, dworkd, n);
  
    for(i=0; i<nrhs; i++){
        j = cublasIzamax( n, dX+i*n, 1) ;
        cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dX+i*n+j-1, 1, &Xnrmv, 1);
        Xnrm = lapackf77_zlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
      
        j = cublasIzamax ( n, dworkd+i*n, 1 );
        cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dworkd+i*n+j-1, 1, &Rnrmv, 1 );
        Rnrm = lapackf77_zlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
      
        if( Rnrm >  Xnrm *cte ){
            goto L10;
        }
    }
    *iter =0; 
    return MAGMA_SUCCESS;
  
  L10:
    ;

    for(iiter=1;iiter<ITERMAX;){
        *info = 0 ; 
        magmablas_zlag2c(n, nrhs, dworkd, n, dSX, n, info );
        if(*info !=0){
            *iter = -2 ;
            goto L40;
        } 
        magma_cpotrs_gpu(uplo, n, nrhs, dSA, ldda, dSX, lddb, info);
      
        for(i=0;i<nrhs;i++){
            magmablas_zcaxpycp(dworks+i*n, dX+i*n, n, dB+i*n,dworkd+i*n) ;
        }
      
        if( nrhs == 1 )
            cublasZhemv(uplo, n, mzone, dA, ldda, dX, 1, zone, dworkd, 1);
        else 
            cublasZhemm(MagmaLeft, uplo, n, nrhs, mzone, dA, ldda, dX, lddx, zone, dworkd, n);
      
        for(i=0; i<nrhs; i++){
            j = cublasIzamax( n, dX+i*n, 1) ;
            cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dX+i*n+j-1, 1, &Xnrmv, 1);
            Xnrm = lapackf77_zlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
          
            j = cublasIzamax ( n, dworkd+i*n, 1 );
            cublasGetMatrix( 1, 1, sizeof(cuDoubleComplex), dworkd+i*n+j-1, 1, &Rnrmv, 1 );
            Rnrm = lapackf77_zlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
          
            if( Rnrm >  Xnrm *cte ){
                goto L20;
            }
        }  

        *iter = iiter;
        return MAGMA_SUCCESS;
      L20:
        iiter++ ;
    }
    *iter = -ITERMAX - 1; 
  
  L40:
    ret = magma_zpotrf_gpu(uplo, n, dA, ldda, info);
    if( (ret != MAGMA_SUCCESS) || (*info != 0) ){
	return ret;
    }
    magmablas_zlacpy( MagmaUpperLower, n, nrhs, dB, lddb, dX, lddx);
    ret = magma_zpotrs_gpu(uplo, n, nrhs, dA, ldda, dX, lddx, info);
    return ret;
}

