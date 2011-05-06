/*
 *  -- MAGMA (version 1.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2010
 *
 * @generated c
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

// Flops formula
#define PRECISION_c
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6. * FMULS_GETRF(m, n) + 2. * FADDS_GETRF(m, n) )
#else
#define FLOPS(m, n) (      FMULS_GETRF(m, n) +      FADDS_GETRF(m, n) )
#endif

float get_LU_error(magma_int_t M, magma_int_t N, 
		    cuFloatComplex *A,  magma_int_t lda, 
		    cuFloatComplex *LU, magma_int_t *IPIV)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    cuFloatComplex alpha = MAGMA_C_ONE;
    cuFloatComplex beta  = MAGMA_C_ZERO;
    cuFloatComplex *L, *U;
    float work[1], matnorm, residual;
                       
    TESTING_MALLOC( L, cuFloatComplex, M*min_mn);
    TESTING_MALLOC( U, cuFloatComplex, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(cuFloatComplex) );
    memset( U, 0, min_mn*N*sizeof(cuFloatComplex) );

    lapackf77_claswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_clacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_clacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_C_MAKE( 1., 0. );
    
    matnorm = lapackf77_clange("f", &M, &N, A, &lda, work);

    blasf77_cgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_C_SUB( LU[i+j*lda], A[i+j*lda] );
	}
    }
    residual = lapackf77_clange("f", &M, &N, LU, &lda, work);

    TESTING_FREE(L);
    TESTING_FREE(U);

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgetrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    TimeStruct       start, end;
    float           flops, gpu_perf, cpu_perf, error;
    cuFloatComplex *h_A, *h_R;
    magma_int_t     *ipiv;

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, lda, ldda;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};

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
        if (M>0 && N>0)
            printf("  testing_cgetrf -M %d -N %d\n\n", M, N);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_cgetrf -M %d -N %d\n\n", 1024, 1024);
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_cgetrf_gpu -M %d -N %d\n\n", 1024, 1024);
        M = N = size[9];
    }
    
    ldda   = ((M+31)/32)*32;
    n2     = M * N;
    min_mn = min(M, N);
    nb     = magma_get_cgetrf_nb(min_mn);

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(ipiv, magma_int_t, min_mn);
    TESTING_MALLOC(    h_A, cuFloatComplex, n2     );
    TESTING_HOSTALLOC( h_R, cuFloatComplex, n2     );

    printf("\n\n");
    printf("  M     N   CPU GFlop/s    GPU GFlop/s   ||PA-LU||/(||A||*N)\n");
    printf("============================================================\n");
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
        lapackf77_clarnv( &ione, ISEED, &n2, h_A );
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_cgetrf(&M, &N, h_A, &lda, ipiv, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of cgetrf had an illegal value.\n", -info);

        cpu_perf = flops / GetTimerValue(start, end);

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_R, &lda, h_A, &lda );
        start = get_current_time();
        magma_cgetrf( M, N, h_R, lda, ipiv, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of cgetrf had an illegal value.\n", -info);

	gpu_perf = flops / GetTimerValue(start, end);

        /* =====================================================================
           Check the factorization
           =================================================================== */
        error = get_LU_error(M, N, h_A, lda, h_R, ipiv);

        printf("%5d %5d  %6.2f         %6.2f         %e\n",
               M, N, cpu_perf, gpu_perf, error);

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( ipiv );
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
