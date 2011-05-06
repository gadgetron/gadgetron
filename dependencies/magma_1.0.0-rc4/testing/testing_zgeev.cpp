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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgeev
*/
#define PRECISION_z

int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    TimeStruct       start, end;
    cuDoubleComplex *h_A, *h_R, *VL, *VR, *h_work, *w1, *w2;
    cuDoubleComplex *w1i, *w2i;
    cuDoubleComplex  mzone = MAGMA_Z_NEG_ONE;
    double          *rwork;
    double           gpu_time, cpu_time, matnorm, result;

    /* Matrix size */
    magma_int_t N=0, n2, lda, nb, lwork;
    magma_int_t size[8] = {1024,2048,3072,4032,5184,6016,7040,8064};

    magma_int_t i, info, checkres, once = 0;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    char *jobl = (char *)"V";
    char *jobr = (char *)"V";

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0) {
                N = atoi(argv[++i]);
                once = 1;
            }
            else if (strcmp("-LN", argv[i])==0)
                jobl = (char *)"N";
            else if (strcmp("-LV", argv[i])==0)
                jobl = (char *)"V";
            else if (strcmp("-RN", argv[i])==0)
                jobr = (char *)"N";
            else if (strcmp("-RV", argv[i])==0)
                jobr = (char *)"V";
        }
        if ( N > 0 )
            printf("  testing_zgeev -L[N|V] -R[N|V] -N %d\n\n", N);
        else
        {
            printf("\nUsage: \n");
            printf("  testing_zgeev -L[N|V] -R[N|V] -N %d\n\n", 1024);
            exit(1);
        }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_zgeev -L[N|V] -R[N|V] -N %d\n\n", 1024);
        N = size[7];
    }

    checkres = getenv("MAGMA_TESTINGS_CHECK") != NULL;

    lda   = N;
    n2    = lda * N;
    nb    = magma_get_zgehrd_nb(N);

#if (defined(PRECISION_s) || defined(PRECISION_d))
    lwork = N*(2+nb);
#else
    lwork = N*(1+nb);
#endif

    w1i   = NULL; 
    w2i   = NULL;
    rwork = NULL;

    TESTING_MALLOC( w1,  cuDoubleComplex, N );
    TESTING_MALLOC( w2,  cuDoubleComplex, N );
#if (defined(PRECISION_s) || defined(PRECISION_d))
    TESTING_MALLOC( w1i, cuDoubleComplex, N );
    TESTING_MALLOC( w2i, cuDoubleComplex, N );
#endif
    TESTING_MALLOC( rwork, double, 2*N ); /* Why is it existing in real ??? */

    TESTING_MALLOC   ( h_A, cuDoubleComplex, n2);
    TESTING_HOSTALLOC( h_R, cuDoubleComplex, n2);
    TESTING_HOSTALLOC( VL , cuDoubleComplex, n2);
    TESTING_HOSTALLOC( VR , cuDoubleComplex, n2);
    TESTING_HOSTALLOC( h_work, cuDoubleComplex, lwork);

    printf("\n\n");
    printf("  N     CPU Time(s)    GPU Time(s)     ||R||_F / ||A||_F\n");
    printf("==========================================================\n");
    for(i=0; i<8; i++){
        if ( argc == 1 ){
            N = size[i];
        }
        
        lda = N;
        n2  = lda*N;

        /* Initialize the matrix */
        lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
#if (defined(PRECISION_c) || defined(PRECISION_z))
        magma_zgeev(jobl[0], jobr[0],
		    N, h_R, lda, w1, 
                    VL, lda, VR, lda,
                    h_work, lwork, rwork, &info);
#else
        magma_zgeev(jobl[0], jobr[0],
 		    N, h_R, lda, w1, w1i,
                    VL, lda, VR, lda,
                    h_work, lwork, &info);
#endif
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of magma_zgeev had an illegal value.\n", -info);

        gpu_time = GetTimerValue(start,end) / 1e3;

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
#if (defined(PRECISION_c) || defined(PRECISION_z))
        lapackf77_zgeev(jobl, jobr,
			&N, h_A, &lda, w2, 
                        VL, &lda, VR, &lda,
			h_work, &lwork, rwork, &info);
#else
        lapackf77_zgeev(jobl, jobr,
			&N, h_A, &lda, w2, w2i, 
                        VL, &lda, VR, &lda,
			h_work, &lwork, &info);
#endif
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of zgeev had an illegal value.\n", -info);

        cpu_time = GetTimerValue(start,end) / 1e3;

        /* =====================================================================
           Check the result compared to LAPACK
           =================================================================== */
        if ( checkres ) {
            matnorm = lapackf77_zlange("f", &N, &ione, w1, &N, rwork);
            printf("norm = %e\n", matnorm);
            blasf77_zaxpy(&N, &mzone, w1, &ione, w2, &ione);

            result = lapackf77_zlange("f", &N, &ione, w2, &N, rwork) / matnorm;

            printf("%5d     %6.2f         %6.2f         %e\n",
                   N, cpu_time, gpu_time, result);
        } else {
            printf("%5d     %6.2f         %6.2f\n",
                   N, cpu_time, gpu_time);
        }

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE(w1);
    TESTING_FREE(w2);
#if (defined(PRECISION_s) || defined(PRECISION_d))
    TESTING_FREE(w1i);
    TESTING_FREE(w2i);
#endif
    TESTING_FREE(rwork);
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );
    TESTING_HOSTFREE( VL  );
    TESTING_HOSTFREE( VR  );
    TESTING_HOSTFREE(rwork);

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
