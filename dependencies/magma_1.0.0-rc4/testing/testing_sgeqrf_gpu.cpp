/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated s

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
#define PRECISION_s
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6.*FMULS_GEQRF(m, n) + 2.*FADDS_GEQRF(m, n) )
#else
#define FLOPS(m, n) (    FMULS_GEQRF(m, n) +    FADDS_GEQRF(m, n) )
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgeqrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    TimeStruct       start, end;
    float           flops, gpu_perf, cpu_perf;
    float           matnorm, work[1];
    float  mzone= MAGMA_S_NEG_ONE;
    float *h_A, *h_R, *tau, *hwork, tmp[1];
    float *d_A;

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, lda, ldda, lhwork;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,9984};

    magma_int_t i, info, min_mn, nb;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0)
                N = atoi(argv[++i]);
            else if (strcmp("-M", argv[i])==0)
                M = atoi(argv[++i]);
        }
        if ( M == 0 ) {
	    M = N;
	}
	if ( N == 0 ) {
	    N = M;
	}
        if (M>0 && N>0)
            printf("  testing_sgeqrf_gpu -M %d -N %d\n\n", M, N);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_sgeqrf_gpu -M %d -N %d\n\n", 1024, 1024);
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_sgeqrf_gpu -M %d -N %d\n\n", 1024, 1024);
        M = N = size[9];
    }

    ldda   = ((M+31)/32)*32;
    n2     = M * N;
    min_mn = min(M, N);
    nb     = magma_get_sgetrf_nb(min_mn);

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(    tau, float, min_mn );
    TESTING_MALLOC(    h_A, float, n2     );
    TESTING_HOSTALLOC( h_R, float, n2     );
    TESTING_DEVALLOC(  d_A, float, ldda*N );

    lhwork = -1;
    lapackf77_sgeqrf(&M, &N, h_A, &M, tau, tmp, &lhwork, &info);
    lhwork = (magma_int_t)MAGMA_S_REAL( tmp[0] );

    TESTING_MALLOC( hwork, float, lhwork );

    printf("\n\n");
    printf("  M     N   CPU GFlop/s   GPU GFlop/s    ||R||_F / ||A||_F\n");
    printf("==========================================================\n");
    for(i=0; i<10; i++){
        if (argc == 1){
	    M = N = size[i];
        }
	min_mn= min(M, N);
	lda   = M;
	n2    = lda*N;
	ldda  = ((M+31)/32)*32;
	flops = FLOPS( (float)M, (float)N ) / 1000000;

        /* Initialize the matrix */
        lapackf77_slarnv( &ione, ISEED, &n2, h_A );
	lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

	/* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_sgeqrf(&M, &N, h_A, &M, tau, hwork, &lhwork, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of lapack_sgeqrf had an illegal value.\n", -info);

        cpu_perf = flops / GetTimerValue(start, end);

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        cublasSetMatrix( M, N, sizeof(float), h_R, lda, d_A, ldda);
        magma_sgeqrf2_gpu( M, N, d_A, ldda, tau, &info);
        cublasSetMatrix( M, N, sizeof(float), h_R, lda, d_A, ldda);

        start = get_current_time();
        magma_sgeqrf2_gpu( M, N, d_A, ldda, tau, &info);
        end = get_current_time();
	if (info < 0)
	  printf("Argument %d of magma_sgeqrf2 had an illegal value.\n", -info);
	
	gpu_perf = flops / GetTimerValue(start, end);
        
        /* =====================================================================
           Check the result compared to LAPACK
           =================================================================== */
        cublasGetMatrix( M, N, sizeof(float), d_A, ldda, h_R, M);
	
        matnorm = lapackf77_slange("f", &M, &N, h_A, &M, work);
        blasf77_saxpy(&n2, &mzone, h_A, &ione, h_R, &ione);
	
        printf("%5d %5d  %6.2f         %6.2f        %e\n",
               M, N, cpu_perf, gpu_perf,
               lapackf77_slange("f", &M, &N, h_R, &M, work) / matnorm);
	
        if (argc != 1)
	  break;
    }
    
    /* Memory clean up */
    TESTING_FREE( tau );
    TESTING_FREE( h_A );
    TESTING_FREE( hwork );
    TESTING_HOSTFREE( h_R );
    TESTING_DEVFREE( d_A );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
