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
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"
#define PRECISION_c

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgesvd
*/
int main( int argc, char** argv) 
{
    cuInit( 0 );
    cublasStatus status = cublasInit();
    if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! CUBLAS initialization error\n");
    }
    printout_devices( );

    cuFloatComplex *h_A, *h_R, *U, *VT, *h_work;
    float *S1, *S2;
#if defined(PRECISION_z) || defined(PRECISION_c)
    float *rwork;
#endif
    float gpu_time, cpu_time;

    TimeStruct start, end;

    /* Matrix size */
    magma_int_t M = 0, N=0, n2, min_mn;
    magma_int_t size[8] = {1024,2048,3072,4032,5184,6016,7040,8064};

    magma_int_t i, j, info;
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
	  printf("  testing_cgesvd -M %d -N %d\n\n", M, N);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_cgesvd -M %d -N %d\n\n", 1024, 1024);

		/* Shutdown */
		status = cublasShutdown();
		if (status != CUBLAS_STATUS_SUCCESS) {
		  fprintf (stderr, "!!!! shutdown error (A)\n");
		}

                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_cgesvd -M %d -N %d\n\n", 1024, 1024);
        N = size[7];
    }

    n2  = M * N;
    min_mn = min(M, N);

    /* Allocate host memory for the matrix */
    h_A = (cuFloatComplex*)malloc(n2 * sizeof(h_A[0]));
    if (h_A == 0) {
        fprintf (stderr, "!!!! host memory allocation error (A)\n");
    }
    VT = (cuFloatComplex*)malloc(N*N * sizeof(h_A[0]));
    if (VT == 0) {
      fprintf (stderr, "!!!! host memory allocation error (VT)\n");
    }
    U = (cuFloatComplex*)malloc(M * M * sizeof(h_A[0]));
    if (U == 0) {
      fprintf (stderr, "!!!! host memory allocation error (U)\n");
    }
    
    S1 = (float*)malloc(min_mn * sizeof(float));
    if (S1 == 0) {
        fprintf (stderr, "!!!! host memory allocation error (S1)\n");
    }
    S2 = (float*)malloc(min_mn * sizeof(float));
    if (S2 == 0) {
      fprintf (stderr, "!!!! host memory allocation error (S2)\n");
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    rwork = (float*)malloc(5 * min_mn * sizeof(float));
    if (rwork == 0) {
      fprintf (stderr, "!!!! host memory allocation error (rwork)\n");
    }
#endif
    cudaMallocHost( (void**)&h_R,  n2*sizeof(cuFloatComplex) );
    if (h_R == 0) {
        fprintf (stderr, "!!!! host memory allocation error (R)\n");
    }

    magma_int_t nb = 128; // magma_get_cgesvd_nb(N);
    magma_int_t lwork = (2*min_mn + max(M,N))*nb;

    cudaMallocHost( (void**)&h_work, lwork*sizeof(cuFloatComplex) );
    if (h_work == 0) {
        fprintf (stderr, "!!!! host memory allocation error (work)\n");
    }

    printf("\n\n");
    printf("  N     CPU Time(s)    GPU Time(s)     ||R||_F / ||A||_F\n");
    printf("==========================================================\n");
    for(i=0; i<8; i++){
        if (argc==1){
            M = N = size[i];
            n2 = M*N;
        }

        /* Initialize the matrix */
        lapackf77_clarnv( &ione, ISEED, &n2, h_A );
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_A, &M, h_R, &M );

#if defined(PRECISION_z) || defined(PRECISION_c)
	magma_cgesvd('A', 'A', M, N,
		     h_R, M, S1, U, M,
		     VT, N, h_work, lwork, rwork, &info); 
#else
	magma_cgesvd('A', 'A', M, N,
		     h_R, M, S1, U, M,
		     VT, N, h_work, lwork, &info); 
#endif
        for(j=0; j<n2; j++)
            h_R[j] = h_A[j];

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
#if defined(PRECISION_z) || defined(PRECISION_c)
	magma_cgesvd('A', 'A', M, N,
                     h_R, M, S1, U, M,
                     VT, N, h_work, lwork, rwork, &info);
#else
	magma_cgesvd('A', 'A', M, N,
		     h_R, M, S1, U, M,
		     VT, N, h_work, lwork, &info); 
#endif
        end = get_current_time();

        gpu_time = GetTimerValue(start,end)/1000.;

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
#if defined(PRECISION_z) || defined(PRECISION_c)
	lapackf77_cgesvd("A", "A", &M, &N,
			 h_A, &M, S2, U, &M,
			 VT, &N, h_work, &lwork, rwork, &info);
#else
	lapackf77_cgesvd("A", "A", &M, &N,
			 h_A, &M, S2, U, &M,
			 VT, &N, h_work, &lwork, &info);
#endif
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of cgesvd had an illegal value.\n", -info);

        cpu_time = GetTimerValue(start,end)/1000.;

        /* =====================================================================
           Check the result compared to LAPACK
           =================================================================== */
        float work[1], matnorm = 1., mone = -1;
        magma_int_t one = 1;

        matnorm = lapackf77_slange("f", &min_mn, &one, S1, &min_mn, work);
        blasf77_saxpy(&min_mn, &mone, S1, &one, S2, &one);

        printf("%5d     %6.2f         %6.2f         %e\n",
               N, cpu_time, gpu_time,
               lapackf77_slange("f", &min_mn, &one, S2, &min_mn, work) / matnorm);

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    free(h_A);
    free(VT);
    free(S1);
    free(S2);
#if defined(PRECISION_z) || defined(PRECISION_c)
    free(rwork);
#endif
    free(U);
    cudaFreeHost(h_work);
    cudaFreeHost(h_R);

    /* Shutdown */
    status = cublasShutdown();
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf (stderr, "!!!! shutdown error (A)\n");
    }
}
