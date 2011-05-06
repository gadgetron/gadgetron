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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgesv_gpu
*/
int main(int argc , char **argv)
{
    TESTING_CUDA_INIT();

    TimeStruct  start, end;
    double      flops, gpu_perf;
    double      Rnorm, Anorm, Bnorm, *work;
    cuDoubleComplex zone  = MAGMA_Z_ONE;
    cuDoubleComplex mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *h_A, *h_B, *h_X;
    cuDoubleComplex *d_A, *d_B;
    magma_int_t *ipiv;
    magma_int_t lda, ldb;
    magma_int_t ldda, lddb;
    magma_int_t i, info, szeA, szeB;
    magma_int_t N        = 0;
    magma_int_t ione     = 1;
    magma_int_t NRHS     = 100;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};
        
    if (argc != 1){
	for(i = 1; i<argc; i++){	
	    if (strcmp("-N", argv[i])==0)
		N = atoi(argv[++i]);
	    else if (strcmp("-nrhs", argv[i])==0)
		NRHS = atoi(argv[++i]);
	}
	if ( N > 0 ) 
	    size[0] = size[9] = N;
    }
    else {
	printf("\nUsage: \n");
	printf("  testing_zgesv_gpu -nrhs %d -N %d\n\n", NRHS, 1024);
    }

    N = size[9];
    ldb = lda = N ;
    lddb = ldda = ((N+31)/32)*32;
    
    TESTING_MALLOC( h_A, cuDoubleComplex, lda*N    );
    TESTING_MALLOC( h_B, cuDoubleComplex, ldb*NRHS );
    TESTING_MALLOC( h_X, cuDoubleComplex, ldb*NRHS );
    TESTING_MALLOC( work, double,         N        );
    TESTING_MALLOC( ipiv, magma_int_t,    N        );

    TESTING_DEVALLOC( d_A, cuDoubleComplex, ldda*N    );
    TESTING_DEVALLOC( d_B, cuDoubleComplex, lddb*NRHS );

    printf("\n\n");
    printf("  N     NRHS       GPU GFlop/s      || b-Ax || / ||A||\n");
    printf("========================================================\n");

    for(i=0; i<10; i++){
	N   = size[i];
        lda = ldb = N;
	ldda = ((N+31)/32)*32;
	lddb = ldda;
	flops = ( FLOPS_GETRF( (double)N, (double)N ) +
		  FLOPS_GETRF( (double)N, (double)NRHS ) ) / 1000000;

        /* Initialize the matrices */
        szeA = lda*N; szeB = ldb*NRHS;
        lapackf77_zlarnv( &ione, ISEED, &szeA, h_A );
        lapackf77_zlarnv( &ione, ISEED, &szeB, h_B );

        cublasSetMatrix( N, N,    sizeof( cuDoubleComplex ), h_A, N, d_A, ldda );
        cublasSetMatrix( N, NRHS, sizeof( cuDoubleComplex ), h_B, N, d_B, lddb );

        //=====================================================================
        // Solve Ax = b through an LU factorization
        //=====================================================================
        start = get_current_time();
        magma_zgesv_gpu( N, NRHS, d_A, ldda, ipiv, d_B, lddb, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of magma_zgesv had an illegal value.\n", -info);

	gpu_perf = flops / GetTimerValue(start, end);

        //=====================================================================
        // ERROR
        //=====================================================================
        cublasGetMatrix( N, NRHS, sizeof( cuDoubleComplex ), d_B, lddb, h_X, ldb );

        Anorm = lapackf77_zlange("I", &N, &N,    h_A, &lda, work);
        Bnorm = lapackf77_zlange("I", &N, &NRHS, h_B, &ldb, work);

        blasf77_zgemm( MagmaNoTransStr, MagmaNoTransStr, &N, &NRHS, &N, 
		       &zone,  h_A, &lda, 
		               h_X, &ldb, 
		       &mzone, h_B, &ldb);
        Rnorm = lapackf77_zlange("I", &N, &NRHS, h_B, &ldb, work);

        printf("%5d  %4d             %6.2f        %e\n",
	       N, NRHS, gpu_perf, Rnorm/(Anorm*Bnorm) );

	if (argc != 1)
	  break;
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( h_B );
    TESTING_FREE( h_X );
    TESTING_FREE( work );
    TESTING_FREE( ipiv );

    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_B );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
