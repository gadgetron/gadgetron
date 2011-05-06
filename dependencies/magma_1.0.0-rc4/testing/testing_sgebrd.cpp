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
#define FLOPS(m, n) ( 6. * FMULS_GEBRD(m, n) + 2. * FADDS_GEBRD(m, n))
#else
#define FLOPS(m, n) (      FMULS_GEBRD(m, n) +      FADDS_GEBRD(m, n))
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgebrd
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    TimeStruct       start, end;
    float           eps, flops, gpu_perf, cpu_perf;
    float *h_A, *h_Q, *h_PT, *h_work, *chkwork;
    float *taup, *tauq;
    float          *diag, *offdiag, *rwork;
    float           result[3] = {0., 0., 0.};

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, lda, lhwork, lchkwork;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};

    magma_int_t i, info, minmn, nb, uselapack, checkres;
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
        if (N>0 && M>0)
            printf("  testing_sgebrd -M %d -N %d\n\n", M, N);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_sgebrd -M %d -N %d\n\n", 1024, 1024);
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_sgebrd -M %d -N %d\n\n", 1024, 1024);
        M = N = size[9];
    }

    uselapack = getenv("MAGMA_USE_LAPACK") != NULL;
    checkres  = getenv("MAGMA_TESTINGS_CHECK") != NULL;

    eps = lapackf77_slamch( "E" );
    lda = M;
    n2  = lda * N;
    nb  = magma_get_sgebrd_nb(N);
    minmn = min(M, N);

    /* Allocate host memory for the matrix */
    TESTING_MALLOC( h_A,     float, lda*N );
    TESTING_MALLOC( tauq,    float, minmn  );
    TESTING_MALLOC( taup,    float, minmn  );
    TESTING_MALLOC( diag,    float, minmn   );
    TESTING_MALLOC( offdiag, float, (minmn-1) );
    TESTING_HOSTALLOC( h_Q, float, lda*N );

    lhwork = (M + N)*nb;
    TESTING_HOSTALLOC( h_work, float, lhwork );

    /* To avoid uninitialized variable warning */
    h_PT    = NULL;
    chkwork = NULL;
    rwork   = NULL; 

    if ( checkres ) {
        lchkwork = max(minmn * nb, M+N);
        /* For optimal performance in sort01 */
        lchkwork = max(lchkwork, minmn*minmn);
        TESTING_MALLOC( h_PT,    float, lda*N   );
        TESTING_MALLOC( chkwork, float, lchkwork );
#if defined(PRECISION_z) || defined(PRECISION_c) 
        TESTING_MALLOC( rwork, float, 5*minmn );
#endif
    }

    printf("\n\n");
    printf("  M    N    CPU GFlop/s    GPU GFlop/s   |A-QHQ'|/N|A|  |I-QQ'|/N \n");
    printf("==================================================================\n");
    for(i=0; i<10; i++){
        if (argc == 1) {
            M = N = size[i];
        }
        minmn = min(M, N);
        lda   = M;
        n2    = lda*N;
        lhwork   = (M + N)*nb;
        lchkwork = max(minmn * nb, M+N);
        /* For optimal performance in sort01 */
        lchkwork = max(lchkwork, minmn*minmn);
        flops = FLOPS( (float)M, (float)N ) / 1e6;

        /* Initialize the matrices */
        lapackf77_slarnv( &ione, ISEED, &n2, h_A );
        lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_Q, &lda );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
        if ( uselapack ) {
            lapackf77_sgebrd( &M, &N, h_Q, &lda, 
                              diag, offdiag, tauq, taup, 
                              h_work, &lhwork, &info);
        } else {
            magma_sgebrd( M, N, h_Q, lda, 
                          diag, offdiag, tauq, taup, 
                          h_work, lhwork, &info);
	}
        end = get_current_time();
        if ( info < 0 )
            printf("Argument %d of lapackf77_sgebrd|magma_sgebrd had an illegal value\n", -info);

        gpu_perf = flops / GetTimerValue(start,end);

        /* =====================================================================
           Check the factorization
           =================================================================== */
        if ( checkres ) {
            lapackf77_slacpy(MagmaUpperLowerStr, &M, &N, h_Q, &lda, h_PT, &lda);
            
            // generate Q & P'
            lapackf77_sorgbr("Q", &M, &minmn, &N, h_Q,  &lda, tauq, chkwork, &lchkwork, &info);
            if ( info < 0 )
              printf("Argument %d of lapackf77_sorgbr had an illegal value\n", -info);
            lapackf77_sorgbr("P", &minmn, &N, &M, h_PT, &lda, taup, chkwork, &lchkwork, &info);
            if ( info < 0 )
              printf("Argument %d of lapackf77_sorgbr (2) had an illegal value\n", -info);
            
            // Test 1:  Check the decomposition A := Q * B * PT
            //      2:  Check the orthogonality of Q
            //      3:  Check the orthogonality of PT
#if defined(PRECISION_z) || defined(PRECISION_c) 
            lapackf77_sbdt01(&M, &N, &ione, 
                             h_A, &lda, h_Q, &lda, 
                             diag, offdiag, h_PT, &lda,
                             chkwork, rwork, &result[0]);
            lapackf77_sort01("Columns", &M, &minmn, h_Q,  &lda, chkwork, &lchkwork, rwork, &result[1]);
            lapackf77_sort01("Rows",    &minmn, &N, h_PT, &lda, chkwork, &lchkwork, rwork, &result[2]);
#else
            lapackf77_sbdt01(&M, &N, &ione, 
                             h_A, &lda, h_Q, &lda, 
                             diag, offdiag, h_PT, &lda,
                             chkwork, &result[0]);
            lapackf77_sort01("Columns", &M, &minmn, h_Q,  &lda, chkwork, &lchkwork, &result[1]);
            lapackf77_sort01("Rows",    &minmn, &N, h_PT, &lda, chkwork, &lchkwork, &result[2]);
#endif
        }

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_sgebrd(&M, &N, h_A, &lda, 
                         diag, offdiag, tauq, taup,
                         h_work, &lhwork, &info);
        end = get_current_time();

        if (info < 0)
            printf("Argument %d of lapackf77_sgebrd had an illegal value.\n", -info);

        cpu_perf = flops / GetTimerValue(start,end);

        /* =====================================================================
           Print performance and error.
           =================================================================== */
        if ( checkres ) {
            printf("%5d %5d   %6.2f        %6.2f       %4.2e %4.2e %4.2e\n",
                   M, N, cpu_perf, gpu_perf,
                   result[0]*eps, result[1]*eps, result[2]*eps );
        } else {
            printf("%5d %5d   %6.2f        %6.2f\n",
                   M, N, cpu_perf, gpu_perf );
        }

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( tauq );
    TESTING_FREE( taup );
    TESTING_FREE( diag );
    TESTING_FREE( offdiag );
    TESTING_HOSTFREE( h_Q );
    TESTING_HOSTFREE( h_work );

    if ( checkres ) {
        TESTING_FREE( h_PT );
        TESTING_FREE( chkwork );
#if defined(PRECISION_z) || defined(PRECISION_c) 
        TESTING_FREE( rwork );
#endif
    }

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
