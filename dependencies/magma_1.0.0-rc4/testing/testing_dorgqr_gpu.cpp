/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d

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
#define PRECISION_d
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n, k) ( 6. * FMULS_UNGQR(m, n, k) + 2. * FADDS_UNGQR(m, n, k))
#else
#define FLOPS(m, n, k) (      FMULS_UNGQR(m, n, k) +      FADDS_UNGQR(m, n, k))
#endif


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dorgqr_gpu
*/
int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    TimeStruct       start, end;
    double           flops, gpu_perf, cpu_perf;
    double           matnorm, work[1];
    double  mzone= MAGMA_D_NEG_ONE;
    double *h_A, *h_R, *tau, *h_work;
    double *d_A, *d_T;

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
            printf("  testing_dorgqr_gpu -M %d -N %d -K %d\n\n", M, N, K);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_dorgqr_gpu  -M %d  -N %d  -K %d\n\n", M, N, K);
                printf("  M, N, and K have to to be K <= N <= M, exit.\n");
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_dorgqr_gpu -M %d  -N %d  -K %d\n\n", 1024, 1024, 1024);
        M = N = K = size[9];
    }
    
    lda    = M;
    ldda   = ((M+31)/32)*32;
    n2     = lda * N;
    min_mn = min(M, N);
    nb     = magma_get_dgeqrf_nb(M);
    lwork  = (M+2*N+nb)*nb;

    TESTING_HOSTALLOC( h_A,    double, lda*N  );
    TESTING_HOSTALLOC( h_work, double, lwork );
    TESTING_MALLOC( h_R, double, lda*N  );
    TESTING_MALLOC( tau, double, min_mn );

    TESTING_DEVALLOC( d_A, double, ldda*N      );
    TESTING_DEVALLOC( d_T, double, ( 2*min_mn+ (N+31)/32*32 )*nb );

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
        nb = magma_get_dgeqrf_nb(M);
        flops = FLOPS( (double)M, (double)N, (double)K ) / 1e6;

        lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
        lapackf77_dlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
        
        cublasSetMatrix( M, N, sizeof(double), h_A, lda, d_A, ldda);
        magma_dgeqrf2_gpu(M, N, d_A, ldda, tau, &info);
        cublasSetMatrix( M, N, sizeof(double), h_A, lda, d_A, ldda);
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magma_dgeqrf_gpu(M, N, d_A, ldda, tau, d_T, &info);
        if ( info < 0)  
            printf("Argument %d of magma_dgeqrf_gpu had an illegal value.\n", -info);
        
        start = get_current_time();
        magma_dorgqr_gpu(M, N, K, d_A, ldda, tau, d_T, nb, &info);
        end = get_current_time();
        if ( info < 0)  
            printf("Argument %d of magma_dorgqr_gpu had an illegal value.\n", -info);
        
        // Get d_A back to the CPU to compare with the CPU result.
        cublasGetMatrix(M, N, sizeof(double), d_A, ldda, h_R, lda);
        
        gpu_perf = flops / GetTimerValue(start,end);
        matnorm = lapackf77_dlange("f", &M, &N, h_A, &lda, work);
        
        /* =====================================================================
           Performs operation using LAPACK 
           =================================================================== */
        lapackf77_dgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
        if ( info < 0)  
            printf("Argument %d of lapackf77_dgeqrf had an illegal value.\n", -info);
        
        start = get_current_time();
        //lapackf77_dorgqr(&M, &N, &K, h_A, &lda, tau, h_work, &lwork, info);
        magma_dorgqr(M, N, K, h_A, lda, tau, d_T, nb, &info);
        end = get_current_time();
        if ( info < 0)  
            printf("Argument %d of magma_dorgqr had an illegal value.\n", -info);
        
        cpu_perf = flops / GetTimerValue(start,end);

        blasf77_daxpy(&n2, &mzone, h_A, &ione, h_R, &ione);
        printf("%5d %5d   %6.1f       %6.1f         %7.2e \n",
               M, N, cpu_perf, gpu_perf,
               lapackf77_dlange("f", &M, &N, h_R, &lda, work) / matnorm );
        
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
