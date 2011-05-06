/*
  -- MAGMA (version 1.0) --
  Univ. of Tennessee, Knoxville
  Univ. of California, Berkeley
  Univ. of Colorado, Denver
  November 2010

  @precisions mixed zc -> ds

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

// Flops formula
#define PRECISION_z
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS_GEQRF(m, n      ) ( 6.*FMULS_GEQRF(m, n      ) + 2.*FADDS_GEQRF(m, n      ) )
#define FLOPS_GEQRS(m, n, nrhs) ( 6.*FMULS_GEQRS(m, n, nrhs) + 2.*FADDS_GEQRS(m, n, nrhs) )
#else
#define FLOPS_GEQRF(m, n      ) (    FMULS_GEQRF(m, n      ) +    FADDS_GEQRF(m, n      ) )
#define FLOPS_GEQRS(m, n, nrhs) (    FMULS_GEQRS(m, n, nrhs) +    FADDS_GEQRS(m, n, nrhs) )
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zcgeqrsv
*/
int main( int argc, char** argv)
{
  /*
#if defined(PRECISION_z) && (GPUSHMEM < 200)
    fprintf(stderr, "This functionnality is not available in MAGMA for this precisions actually\n");
    return EXIT_SUCCESS;
#else
  */
    TESTING_CUDA_INIT();

    TimeStruct  start, end;
    double      flops, gpu_perf, cpu_perf;
    double      gpu_perfd, gpu_perfs;
    double      Rnorm, Anorm, work[1];
    cuDoubleComplex zone  = MAGMA_Z_ONE;
    cuDoubleComplex mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *h_A, *h_B, *h_X, *h_R;
    cuDoubleComplex *d_A, *d_B, *d_X, *d_T;
    cuFloatComplex  *d_SA, *d_SB, *d_ST;
    cuDoubleComplex *h_workd, *tau, tmp[1];
    cuFloatComplex  *h_works, *tau_s;
    magma_int_t lda,  ldb,  ldx, lhwork;
    magma_int_t ldda, lddb, lddx;
    magma_int_t i, iter, info, size, min_mn, nb;
    magma_int_t M        = 0;
    magma_int_t N        = 0;
    magma_int_t ione     = 1;
    magma_int_t NRHS     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t sizetest[10] = {1024,2048,3072,4032,5184,6016,7040,7520,8064,8192};
        
    if (argc != 1){
	for(i = 1; i<argc; i++){	
	    if (strcmp("-N", argv[i])==0)
		N = atoi(argv[++i]);
	    else if (strcmp("-M", argv[i])==0)
		M = atoi(argv[++i]);
	    else if (strcmp("-nrhs", argv[i])==0)
		NRHS = atoi(argv[++i]);
	}
	if (N>0) sizetest[0] = sizetest[9] = N;
	else exit(1);
    }
    else {
	printf("\nUsage: \n");
	printf("  testing_zcgeqrsv_gpu -M %d -N %d -nrhs %d\n\n", 1024, 1024, NRHS);
    }
    printf("Epsilon(double): %8.6e\n"
	   "Epsilon(single): %8.6e\n\n", 
	   lapackf77_dlamch("Epsilon"), lapackf77_slamch("Epsilon") );

    N = sizetest[9];
    if (M == 0) M = N;
    ldb  = lda = M;
    ldx  = N;
    ldda = ((M+31)/32)*32;
    lddb = ldda;
    lddx = ((N+31)/32)*32;    
    nb = max( magma_get_zgeqrf_nb( M ), magma_get_cgeqrf_nb( M ) );
    min_mn = min( M, N );

    TESTING_MALLOC( h_A, cuDoubleComplex, lda*N    );
    TESTING_MALLOC( h_B, cuDoubleComplex, ldb*NRHS );
    TESTING_MALLOC( h_X, cuDoubleComplex, ldx*NRHS );
    TESTING_MALLOC( h_R, cuDoubleComplex, ldb*NRHS );
    TESTING_MALLOC( tau, cuDoubleComplex, min_mn   );

    TESTING_DEVALLOC( d_A, cuDoubleComplex, ldda*N      );
    TESTING_DEVALLOC( d_B, cuDoubleComplex, lddb*NRHS   );
    TESTING_DEVALLOC( d_X, cuDoubleComplex, lddx*NRHS   );
    TESTING_DEVALLOC( d_T, cuDoubleComplex, ( 2*min_mn+ (N+31)/32*32 )*nb );
    tau_s = (cuFloatComplex*)tau;
    d_ST  = (cuFloatComplex*)d_T;

    lhwork = -1;
    lapackf77_zgels( MagmaNoTransStr, &M, &N, &NRHS,
		     h_A, &lda, h_B, &ldb, tmp, &lhwork, &info);
    lhwork = (magma_int_t)MAGMA_Z_REAL( tmp[0] );
    lhwork = max( lhwork, (nb * max((M-N+nb+2*(NRHS)), 1) ) );

    TESTING_MALLOC( h_workd, cuDoubleComplex, lhwork );
    h_works = (cuFloatComplex*)h_workd;

    printf("\n\n");
    printf("        CPU GFlop/s         G P U  GFlop/s   \n");
    printf("  N         DP          DP       SP       MP    ||b-Ax||/||A||  NumIter\n");
    printf("=======================================================================\n");

    for(i=0; i<8; i++){
        M = N = sizetest[i];
	lda = ldb = M;
	ldx = N;
	ldda = lddb = ((M+31)/32) * 32;
	lddx = ((N+31)/32) * 32;

	flops = ( FLOPS_GEQRF( (double)M, (double)N ) 
		  + FLOPS_GEQRS( (double)M, (double)N, (double)NRHS ) ) / 1000000;

	/* Initialize matrices */
        size = lda*N;
        lapackf77_zlarnv( &ione, ISEED, &size, h_A );
        size = ldb*NRHS;
        lapackf77_zlarnv( &ione, ISEED, &size, h_B );
        lapackf77_zlacpy( MagmaUpperLowerStr, &M, &NRHS, h_B, &ldb, h_R, &ldb );

        cublasSetMatrix( M, N,    sizeof(cuDoubleComplex), h_A, lda, d_A, ldda );
        cublasSetMatrix( M, NRHS, sizeof(cuDoubleComplex), h_B, ldb, d_B, lddb );

        //=====================================================================
        //              Mixed Precision Iterative Refinement - GPU
        //=====================================================================
        start = get_current_time();
        magma_zcgeqrsv_gpu( M, N, NRHS, 
                            d_A, ldda, d_B, lddb, 
			    d_X, lddx, &iter, &info );
        end = get_current_time();
	if (info < 0)
	  printf("Argument %d of magma_zcgeqrsv had an illegal value.\n", -info);
	gpu_perf = flops / GetTimerValue(start, end);
        
        //=====================================================================
        //                 Error Computation
        //=====================================================================
        cublasGetMatrix(N, NRHS, sizeof(cuDoubleComplex), d_X, lddx, h_X, ldx);
        blasf77_zgemm( MagmaNoTransStr, MagmaNoTransStr,
		       &M, &NRHS, &N,
                       &mzone, h_A, &lda, 
                               h_X, &ldx, 
                       &zone,  h_R, &ldb);
        Anorm = lapackf77_zlange("f", &N, &N,    h_A, &lda, work);
        Rnorm = lapackf77_zlange("f", &N, &NRHS, h_R, &ldb, work);

        //=====================================================================
        //                 Double Precision Solve
        //=====================================================================
        cublasSetMatrix( M, N,    sizeof(cuDoubleComplex), h_A, lda, d_A, ldda );
        cublasSetMatrix( M, NRHS, sizeof(cuDoubleComplex), h_B, ldb, d_B, lddb );

	start = get_current_time();
        magma_zgels_gpu( MagmaNoTrans, M, N, NRHS, d_A, ldda, 
			 d_B, lddb, h_workd, lhwork, &info);
        end = get_current_time();
	gpu_perfd = flops / GetTimerValue(start, end);


        //=====================================================================
        //                 Single Precision Solve
        //=====================================================================
        cublasSetMatrix( M, N,    sizeof(cuDoubleComplex), h_A, lda, d_A, ldda );
        cublasSetMatrix( M, NRHS, sizeof(cuDoubleComplex), h_B, ldb, d_B, lddb );

	/* The allocation of d_SA and d_SB is done here to avoid 
	   to double the memory used on GPU with zcgeqrsv */
	TESTING_DEVALLOC( d_SA, cuFloatComplex, ldda*N      );
	TESTING_DEVALLOC( d_SB, cuFloatComplex, lddb*NRHS   );
	magmablas_zlag2c(M, N,    d_A, ldda, d_SA, ldda, &info);
	magmablas_zlag2c(N, NRHS, d_B, lddb, d_SB, lddb, &info);

	start = get_current_time();
        magma_cgels_gpu( MagmaNoTrans, M, N, NRHS, d_SA, ldda, 
			 d_SB, lddb, h_works, lhwork, &info);
        end = get_current_time();
	gpu_perfs = flops / GetTimerValue(start, end);
	TESTING_DEVFREE( d_SA );
	TESTING_DEVFREE( d_SB );

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_zgels( MagmaNoTransStr, &M, &N, &NRHS,
			 h_A, &lda, h_B, &ldb, h_workd, &lhwork, &info);
        end = get_current_time();
        cpu_perf = flops / GetTimerValue(start, end);

        printf("%5d  %8.2f   %9.2f   %6.2f   %6.2f    %e    %2d\n",
               sizetest[i], cpu_perf, gpu_perfd, gpu_perfs, gpu_perf, Rnorm / Anorm, iter );

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( h_B );
    TESTING_FREE( h_X );
    TESTING_FREE( h_R );
    TESTING_FREE( tau );
    TESTING_FREE( h_workd );

    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_B );
    TESTING_DEVFREE( d_X );
    TESTING_DEVFREE( d_T );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();

    //#endif /*defined(PRECISION_z) && (GPUSHMEM < 200)*/
}
