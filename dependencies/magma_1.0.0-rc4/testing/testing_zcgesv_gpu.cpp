/*
  -- MAGMA (version 1.0) --
  Univ. of Tennessee, Knoxville
  Univ. of California, Berkeley
  Univ. of Colorado, Denver
  November 2010

  @precisions mixed zc -> ds

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define PRECISION_z
// Flops formula
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS_GETRF(m, n   ) ( 6.*FMULS_GETRF(m, n   ) + 2.*FADDS_GETRF(m, n   ) )
#define FLOPS_GETRS(m, nrhs) ( 6.*FMULS_GETRS(m, nrhs) + 2.*FADDS_GETRS(m, nrhs) )
#else
#define FLOPS_GETRF(m, n   ) (    FMULS_GETRF(m, n   ) +    FADDS_GETRF(m, n   ) )
#define FLOPS_GETRS(m, nrhs) (    FMULS_GETRS(m, nrhs) +    FADDS_GETRS(m, nrhs) )
#endif

int main(int argc , char **argv)
{
    TESTING_CUDA_INIT();

    TimeStruct  start, end;
    double      flopsF, flopsS, gpu_perf;
    double      gpu_perfdf, gpu_perfds;
    double      gpu_perfsf, gpu_perfss;
    double      Rnorm, Anorm;
    cuDoubleComplex zone  = MAGMA_Z_ONE;
    cuDoubleComplex mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *h_A, *h_B, *h_X;
    cuDoubleComplex *d_A, *d_B, *d_X, *d_WORKD;
    cuFloatComplex  *d_As, *d_Bs, *d_WORKS;
    double          *h_workd;
    magma_int_t *h_ipiv, *d_ipiv;
    magma_int_t lda, ldb, ldx;
    magma_int_t ldda, lddb, lddx;
    magma_int_t i, iter, info, size;
    magma_int_t N        = 0;
    magma_int_t ione     = 1;
    magma_int_t NRHS     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t sizetest[10] = {1024,2048,3072,4032,5184,6016,7040,7520,8064,8192};

    char trans = MagmaNoTrans;
    //char trans = MagmaConjTrans;
    char trans_str[2] = {trans, 0};

    if (argc != 1){
	for(i = 1; i<argc; i++){	
	    if (strcmp("-N", argv[i])==0)
		N = atoi(argv[++i]);
	    else if (strcmp("-nrhs", argv[i])==0)
		NRHS = atoi(argv[++i]);
	}
	if (N>0) sizetest[0] = sizetest[9] = N;
	else exit(1);
    }
    else {
	printf("\nUsage: \n");
	printf("  testing_zcgesv_gpu -nrhs %d -N %d\n\n", NRHS, 1024);
    }
    printf("Epsilon(double): %8.6e\n"
	   "Epsilon(single): %8.6e\n\n", 
	   lapackf77_dlamch("Epsilon"), lapackf77_slamch("Epsilon") );

    N = sizetest[9];
    ldb  = ldx = lda = N;
    ldda = ((N+31)/32)*32;
    lddb = lddx = ldda;

    TESTING_MALLOC( h_A, cuDoubleComplex, lda*N    );
    TESTING_MALLOC( h_B, cuDoubleComplex, ldb*NRHS );
    TESTING_MALLOC( h_X, cuDoubleComplex, ldx*NRHS );
    TESTING_MALLOC( h_ipiv, magma_int_t,    N        );
    TESTING_MALLOC( h_workd, double, N );
    
    TESTING_DEVALLOC( d_A,     cuDoubleComplex, ldda*N        );
    TESTING_DEVALLOC( d_B,     cuDoubleComplex, lddb*NRHS     );
    TESTING_DEVALLOC( d_X,     cuDoubleComplex, lddx*NRHS     );
    TESTING_DEVALLOC( d_ipiv,  magma_int_t,     N             );
    TESTING_DEVALLOC( d_WORKS, cuFloatComplex,  ldda*(N+NRHS) );
    TESTING_DEVALLOC( d_WORKD, cuDoubleComplex, N*NRHS        );

    printf("  N   DP-Factor  DP-Solve  SP-Factor  SP-Solve  MP-Solve  ||b-Ax||/||A||  NumIter\n");
    printf("==================================================================================\n");
    for(i=0; i<10; i++){
	N = sizetest[i] ;
	
	flopsF = FLOPS_GETRF( (double)N, (double)N ) / 1000000;
	flopsS = flopsF + ( FLOPS_GETRS( (double)N, (double)NRHS ) / 1000000 );

	ldb  = ldx = lda = N;
	ldda = ((N+31)/32)*32;
	lddb = lddx = N;//ldda;

	/* Initialize matrices */
        size = lda * N;
        lapackf77_zlarnv( &ione, ISEED, &size, h_A );
        size = ldb * NRHS;
        lapackf77_zlarnv( &ione, ISEED, &size, h_B );
        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &NRHS, h_B, &ldb, h_X, &ldx);

        printf("%5d  ",N);

        cublasSetMatrix( N, N,    sizeof( cuDoubleComplex ), h_A, lda, d_A, ldda );
        cublasSetMatrix( N, NRHS, sizeof( cuDoubleComplex ), h_B, ldb, d_B, lddb );

        //=====================================================================
        //              MIXED - GPU
        //=====================================================================
        start = get_current_time();
        magma_zcgesv_gpu( trans, N, NRHS, 
			  d_A, ldda, h_ipiv, d_ipiv, 
			  d_B, lddb, d_X, lddx, 
			  d_WORKD, d_WORKS, &iter, &info);
	end = get_current_time();
	if (info < 0)
            printf("Argument %d of magma_zcgesv had an illegal value.\n", -info);
	gpu_perf = flopsS / GetTimerValue(start, end);

        //=====================================================================
        //              ERROR DP vs MIXED  - GPU
        //=====================================================================
        cublasGetMatrix( N, NRHS, sizeof( cuDoubleComplex ), d_X, lddx, h_X, ldx );

        Anorm = lapackf77_zlange("I", &N, &N, h_A, &lda, h_workd);
        blasf77_zgemm( trans_str, MagmaNoTransStr, 
                       &N, &NRHS, &N, 
                       &zone,  h_A, &lda,
                               h_X, &ldx,
                       &mzone, h_B, &ldb);
        Rnorm = lapackf77_zlange("I", &N, &NRHS, h_B, &ldb, h_workd);

        //=====================================================================
	//                 Double Precision Factor 
	//=====================================================================
	cublasSetMatrix( N, N, sizeof(cuDoubleComplex), h_A, lda, d_A, ldda );
        
        start = get_current_time();
        magma_zgetrf_gpu(N, N, d_A, ldda, h_ipiv, &info);
        end = get_current_time();
	if (info < 0)
	    printf("Argument %d of magma_zgetrf had an illegal value.\n", -info);
	gpu_perfdf = flopsF / GetTimerValue(start, end);

	printf("%6.2f    ", gpu_perfdf); fflush(stdout);

        //=====================================================================
	//                 Double Precision Solve 
	//=====================================================================
	cublasSetMatrix( N, N,    sizeof(cuDoubleComplex), h_A, lda, d_A, ldda );
	cublasSetMatrix( N, NRHS, sizeof(cuDoubleComplex), h_B, ldb, d_B, lddb );

        start = get_current_time();
        magma_zgetrf_gpu(N, N, d_A, ldda, h_ipiv, &info);
        magma_zgetrs_gpu( trans, N, NRHS, d_A, ldda, h_ipiv, d_B, lddb, &info );
        end = get_current_time();
	if (info < 0)
	    printf("Argument %d of magma_zgetrs had an illegal value.\n", -info);

	gpu_perfds = flopsS / GetTimerValue(start, end);

	printf("%6.2f    ", gpu_perfds); fflush(stdout);

        //=====================================================================
	//                 Single Precision Factor 
	//=====================================================================
	d_As = d_WORKS;
	d_Bs = d_WORKS + ldda*N;
	cublasSetMatrix( N, N,    sizeof(cuDoubleComplex), h_A, lda, d_A, ldda );
	cublasSetMatrix( N, NRHS, sizeof(cuDoubleComplex), h_B, ldb, d_B, lddb );
	magmablas_zlag2c(N, N,    d_A, ldda, d_As, ldda, &info ); 
	magmablas_zlag2c(N, NRHS, d_B, lddb, d_Bs, lddb, &info );

        start = get_current_time();
        magma_cgetrf_gpu(N, N, d_As, ldda, h_ipiv, &info);
        end = get_current_time();
	if (info < 0)
	    printf("Argument %d of magma_cgetrf had an illegal value.\n", -info);
	
	gpu_perfsf = flopsF / GetTimerValue(start, end);
	printf("%6.2f     ", gpu_perfsf); fflush(stdout);

	//=====================================================================
	//                 Single Precision Solve 
	//=====================================================================
	magmablas_zlag2c(N, N,    d_A, ldda, d_As, ldda, &info ); 
	magmablas_zlag2c(N, NRHS, d_B, lddb, d_Bs, lddb, &info );

        start = get_current_time();
        magma_cgetrf_gpu( N, N,    d_As, ldda, h_ipiv, &info);
        magma_cgetrs_gpu( trans, N, NRHS, d_As, ldda, h_ipiv, 
			  d_Bs, lddb, &info);
        end = get_current_time();
        if (info < 0)
	    printf("Argument %d of magma_cgetrs had an illegal value.\n", -info);

	gpu_perfss = flopsS / GetTimerValue(start, end);
	printf("%6.2f     ", gpu_perfss); fflush(stdout);

	printf("%6.2f     ", gpu_perf);
	printf("%e    %3d\n", Rnorm/Anorm, iter); fflush(stdout);

        if( argc != 1 ) break ;        
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( h_B );
    TESTING_FREE( h_X );
    TESTING_FREE( h_ipiv );
    TESTING_FREE( h_workd );
    
    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_B );
    TESTING_DEVFREE( d_X );
    TESTING_DEVFREE( d_ipiv );
    TESTING_DEVFREE( d_WORKS );
    TESTING_DEVFREE( d_WORKD );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
