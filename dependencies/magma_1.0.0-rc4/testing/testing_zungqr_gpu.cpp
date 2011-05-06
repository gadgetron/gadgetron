/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

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
#define FLOPS(m, n, k) ( 6. * FMULS_UNGQR(m, n, k) + 2. * FADDS_UNGQR(m, n, k))
#else
#define FLOPS(m, n, k) (      FMULS_UNGQR(m, n, k) +      FADDS_UNGQR(m, n, k))
#endif


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zungqr_gpu
*/
int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    TimeStruct       start, end;
    double           flops, gpu_perf, cpu_perf;
    double           matnorm, work[1];
    cuDoubleComplex  mzone= MAGMA_Z_NEG_ONE;
    cuDoubleComplex *h_A, *h_R, *tau, *h_work;
    cuDoubleComplex *d_A, *d_T;

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
            printf("  testing_zungqr_gpu -M %d -N %d -K %d\n\n", M, N, K);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_zungqr_gpu  -M %d  -N %d  -K %d\n\n", M, N, K);
                printf("  M, N, and K have to to be K <= N <= M, exit.\n");
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_zungqr_gpu -M %d  -N %d  -K %d\n\n", 1024, 1024, 1024);
        M = N = K = size[9];
    }
    
    lda    = M;
    ldda   = ((M+31)/32)*32;
    n2     = lda * N;
    min_mn = min(M, N);
    nb     = magma_get_zgeqrf_nb(M);
    lwork  = (M+2*N+nb)*nb;

    TESTING_HOSTALLOC( h_A,    cuDoubleComplex, lda*N  );
    TESTING_HOSTALLOC( h_work, cuDoubleComplex, lwork );
    TESTING_MALLOC( h_R, cuDoubleComplex, lda*N  );
    TESTING_MALLOC( tau, cuDoubleComplex, min_mn );

    TESTING_DEVALLOC( d_A, cuDoubleComplex, ldda*N      );
    TESTING_DEVALLOC( d_T, cuDoubleComplex, ( 2*min_mn+ (N+31)/32*32 )*nb );

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
        nb = magma_get_zgeqrf_nb(M);
        flops = FLOPS( (double)M, (double)N, (double)K ) / 1e6;

        lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
        lapackf77_zlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
        
        cublasSetMatrix( M, N, sizeof(cuDoubleComplex), h_A, lda, d_A, ldda);
        magma_zgeqrf2_gpu(M, N, d_A, ldda, tau, &info);
        cublasSetMatrix( M, N, sizeof(cuDoubleComplex), h_A, lda, d_A, ldda);
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magma_zgeqrf_gpu(M, N, d_A, ldda, tau, d_T, &info);
        if ( info < 0)  
            printf("Argument %d of magma_zgeqrf_gpu had an illegal value.\n", -info);
        
        start = get_current_time();
        magma_zungqr_gpu(M, N, K, d_A, ldda, tau, d_T, nb, &info);
        end = get_current_time();
        if ( info < 0)  
            printf("Argument %d of magma_zungqr_gpu had an illegal value.\n", -info);
        
        // Get d_A back to the CPU to compare with the CPU result.
        cublasGetMatrix(M, N, sizeof(cuDoubleComplex), d_A, ldda, h_R, lda);
        
        gpu_perf = flops / GetTimerValue(start,end);
        matnorm = lapackf77_zlange("f", &M, &N, h_A, &lda, work);
        
        /* =====================================================================
           Performs operation using LAPACK 
           =================================================================== */
        lapackf77_zgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
        if ( info < 0)  
            printf("Argument %d of lapackf77_zgeqrf had an illegal value.\n", -info);
        
        start = get_current_time();
        //lapackf77_zungqr(&M, &N, &K, h_A, &lda, tau, h_work, &lwork, info);
        magma_zungqr(M, N, K, h_A, lda, tau, d_T, nb, &info);
        end = get_current_time();
        if ( info < 0)  
            printf("Argument %d of magma_zungqr had an illegal value.\n", -info);
        
        cpu_perf = flops / GetTimerValue(start,end);

        blasf77_zaxpy(&n2, &mzone, h_A, &ione, h_R, &ione);
        printf("%5d %5d   %6.1f       %6.1f         %7.2e \n",
               M, N, cpu_perf, gpu_perf,
               lapackf77_zlange("f", &M, &N, h_R, &lda, work) / matnorm );
        
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
