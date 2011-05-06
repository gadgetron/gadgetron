/*
 *  -- MAGMA (version 1.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2010
 *
 * @precisions normal z -> c d s
 *
 **/
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

#define PRECISION_z
// Flops formula
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_POTRF(n) + 2. * FADDS_POTRF(n) )
#else
#define FLOPS(n) (      FMULS_POTRF(n) +      FADDS_POTRF(n) )
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zpotrf
*/
int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    TimeStruct  start, end;
    double      flops, gpu_perf, cpu_perf;
    cuDoubleComplex *h_A, *h_R;
    cuDoubleComplex *d_A;
    magma_int_t N = 0, n2, lda, ldda;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6048,7200,8064,8928,10240};
    
    magma_int_t i, info;
    const char *uplo     = MagmaUpperStr;
    cuDoubleComplex mzone= MAGMA_Z_NEG_ONE;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    double      work[1], matnorm;
    
    if (argc != 1){
	for(i = 1; i<argc; i++){	
	    if (strcmp("-N", argv[i])==0)
		N = atoi(argv[++i]);
	}
	if (N>0) size[0] = size[9] = N;
	else exit(1);
    }
    else {
	printf("\nUsage: \n");
	printf("  testing_zpotrf_gpu -N %d\n\n", 1024);
    }

    /* Allocate host memory for the matrix */
    n2   = size[9] * size[9];
    ldda = ((size[9]+31)/32) * 32;
    TESTING_MALLOC(    h_A, cuDoubleComplex, n2);
    TESTING_HOSTALLOC( h_R, cuDoubleComplex, n2);
    TESTING_DEVALLOC(  d_A, cuDoubleComplex, ldda*size[9] );

    printf("\n\n");
    printf("  N    CPU GFlop/s    GPU GFlop/s    ||R||_F / ||A||_F\n");
    printf("========================================================\n");
    for(i=0; i<10; i++){
	N   = size[i];
	lda = N; 
	n2  = lda*N;
        flops = FLOPS( (double)N ) / 1000000;
	
	ldda = ((N+31)/32)*32;

	/* Initialize the matrix */
	lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
        /* Symmetrize and increase the diagonal */
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_Z_SET2REAL( h_A[i*lda+i], ( MAGMA_Z_GET_X(h_A[i*lda+i]) + 1.*N ) );
                for(j=0; j<i; j++)
                    h_A[i*lda+j] = cuConj(h_A[j*lda+i]);
            }
        }
        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );

	/* ====================================================================
	   Performs operation using MAGMA 
	   =================================================================== */
	cublasSetMatrix( N, N, sizeof(cuDoubleComplex), h_A, lda, d_A, ldda);
	magma_zpotrf_gpu(uplo[0], N, d_A, ldda, &info);

	cublasSetMatrix( N, N, sizeof(cuDoubleComplex), h_A, lda, d_A, ldda);
      	start = get_current_time();
	magma_zpotrf_gpu(uplo[0], N, d_A, ldda, &info);
	end = get_current_time();
	if (info < 0)
            printf("Argument %d of magma_zpotrf had an illegal value.\n", -info);

        gpu_perf = flops / GetTimerValue(start, end);
	
	/* =====================================================================
	   Performs operation using LAPACK 
	   =================================================================== */
	start = get_current_time();
	lapackf77_zpotrf(uplo, &N, h_A, &lda, &info);
	end = get_current_time();
	if (info < 0)  
	    printf("Argument %d of zpotrf had an illegal value.\n", -info);
	
        cpu_perf = flops / GetTimerValue(start, end);
      
	/* =====================================================================
	   Check the result compared to LAPACK
	   =================================================================== */
	cublasGetMatrix( N, N, sizeof(cuDoubleComplex), d_A, ldda, h_R, lda);
	matnorm = lapackf77_zlange("f", &N, &N, h_A, &lda, work);
	blasf77_zaxpy(&n2, &mzone, h_A, &ione, h_R, &ione);
	printf("%5d    %6.2f         %6.2f        %e\n", 
	       size[i], cpu_perf, gpu_perf,
	       lapackf77_zlange("f", &N, &N, h_R, &lda, work) / matnorm);
	
	if (argc != 1)
	    break;
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );
    TESTING_DEVFREE( d_A );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
