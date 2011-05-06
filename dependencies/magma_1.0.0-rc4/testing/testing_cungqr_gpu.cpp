/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated c

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
#define PRECISION_c
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n, k) ( 6. * FMULS_UNGQR(m, n, k) + 2. * FADDS_UNGQR(m, n, k))
#else
#define FLOPS(m, n, k) (      FMULS_UNGQR(m, n, k) +      FADDS_UNGQR(m, n, k))
#endif


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cungqr_gpu
*/
int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    TimeStruct       start, end;
    float           flops, gpu_perf, cpu_perf;
    float           matnorm, work[1];
    cuFloatComplex  mzone= MAGMA_C_NEG_ONE;
    cuFloatComplex *h_A, *h_R, *tau, *h_work;
    cuFloatComplex *d_A, *d_T;

    /* Matrix size */
    magma_int_t M=0, N=0, K=0;
    magma_int_t n2, lda, ldda, lwork, min_mn, nb;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,9984};
    
    magma_int_t i, info;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    if (argc != 1){
        for(i = 1; i<argc; i++){	
            if (strcmp("-N", argv[i])==0)
                N = atoi(argv[++i]);
            else if (strcmp("-M", argv[i])==0)
                M = atoi(argv[++i]);
            else if (strcmp("-K", argv[i])==0)
                K = atoi(argv[++i]);
        }
        if (N>0 && M>0 && M >= N && K >0 && K <= N)
            printf("  testing_cungqr_gpu -M %d -N %d -K %d\n\n", M, N, K);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_cungqr_gpu  -M %d  -N %d  -K %d\n\n", M, N, K);
                printf("  M, N, and K have to to be K <= N <= M, exit.\n");
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_cungqr_gpu -M %d  -N %d  -K %d\n\n", 1024, 1024, 1024);
        M = N = K = size[9];
    }
    
    lda    = M;
    ldda   = ((M+31)/32)*32;
    n2     = lda * N;
    min_mn = min(M, N);
    nb     = magma_get_cgeqrf_nb(M);
    lwork  = (M+2*N+nb)*nb;

    TESTING_HOSTALLOC( h_A,    cuFloatComplex, lda*N  );
    TESTING_HOSTALLOC( h_work, cuFloatComplex, lwork );
    TESTING_MALLOC( h_R, cuFloatComplex, lda*N  );
    TESTING_MALLOC( tau, cuFloatComplex, min_mn );

    TESTING_DEVALLOC( d_A, cuFloatComplex, ldda*N      );
    TESTING_DEVALLOC( d_T, cuFloatComplex, ( 2*min_mn+ (N+31)/32*32 )*nb );

    printf("\n");
    printf("  M     N    CPU GFlop/s   GPU GFlop/s   ||R|| / ||A||\n");
    printf("=======================================================\n");
    for(i=0; i<10; i++){
        if (argc == 1){
            M = N = size[i];
            K = min(M, N);
        }
        lda  = M;
        ldda = ((M+31)/32)*32;
        n2 = lda*N;
        nb = magma_get_cgeqrf_nb(M);
        flops = FLOPS( (float)M, (float)N, (float)K ) / 1e6;

        lapackf77_clarnv( &ione, ISEED, &n2, h_A );
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
        
        cublasSetMatrix( M, N, sizeof(cuFloatComplex), h_A, lda, d_A, ldda);
        magma_cgeqrf2_gpu(M, N, d_A, ldda, tau, &info);
        cublasSetMatrix( M, N, sizeof(cuFloatComplex), h_A, lda, d_A, ldda);
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magma_cgeqrf_gpu(M, N, d_A, ldda, tau, d_T, &info);
        if ( info < 0)  
            printf("Argument %d of magma_cgeqrf_gpu had an illegal value.\n", -info);
        
        start = get_current_time();
        magma_cungqr_gpu(M, N, K, d_A, ldda, tau, d_T, nb, &info);
        end = get_current_time();
        if ( info < 0)  
            printf("Argument %d of magma_cungqr_gpu had an illegal value.\n", -info);
        
        // Get d_A back to the CPU to compare with the CPU result.
        cublasGetMatrix(M, N, sizeof(cuFloatComplex), d_A, ldda, h_R, lda);
        
        gpu_perf = flops / GetTimerValue(start,end);
        matnorm = lapackf77_clange("f", &M, &N, h_A, &lda, work);
        
        /* =====================================================================
           Performs operation using LAPACK 
           =================================================================== */
        lapackf77_cgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
        if ( info < 0)  
            printf("Argument %d of lapackf77_cgeqrf had an illegal value.\n", -info);
        
        start = get_current_time();
        //lapackf77_cungqr(&M, &N, &K, h_A, &lda, tau, h_work, &lwork, info);
        magma_cungqr(M, N, K, h_A, lda, tau, d_T, nb, &info);
        end = get_current_time();
        if ( info < 0)  
            printf("Argument %d of magma_cungqr had an illegal value.\n", -info);
        
        cpu_perf = flops / GetTimerValue(start,end);

        blasf77_caxpy(&n2, &mzone, h_A, &ione, h_R, &ione);
        printf("%5d %5d   %6.1f       %6.1f         %7.2e \n",
               M, N, cpu_perf, gpu_perf,
               lapackf77_clange("f", &M, &N, h_R, &lda, work) / matnorm );
        
        if (argc != 1)
            break;
    }
    
    /* Memory clean up */
    TESTING_HOSTFREE( h_A );
    TESTING_HOSTFREE( h_work );
    TESTING_FREE( h_R );
    TESTING_FREE( tau );

    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_T );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
