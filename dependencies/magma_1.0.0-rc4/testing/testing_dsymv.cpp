/*
 *  -- MAGMA (version 1.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2010
 *
 *  @generated d
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cblas.h>

#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define PRECISION_d
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_SYMV(n) + 2. * FADDS_SYMV(n))
#else
#define FLOPS(n) (      FMULS_SYMV(n) +      FADDS_SYMV(n))
#endif

int main(int argc, char **argv)
{	
    TESTING_CUDA_INIT();

    TimeStruct  start, end;
    double      flops, magma_perf, cuda_perf, error, work[1];
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    double mzone = MAGMA_D_NEG_ONE;

    FILE        *fp ; 
    magma_int_t N, m, i, lda, LDA;
    magma_int_t matsize;
    magma_int_t vecsize;
    magma_int_t istart = 64;
    magma_int_t incx = 1;
    char        uplo = MagmaLower;
    double alpha = MAGMA_D_MAKE(1., 0.); // MAGMA_D_MAKE(  1.5, -2.3 );
    double beta  = MAGMA_D_MAKE(0., 0.); // MAGMA_D_MAKE( -0.6,  0.8 );
    double *A, *X, *Y, *Ycublas, *Ymagma;
    double *dA, *dX, *dY;
#if defined(PRECISION_z) || defined(PRECISION_c)

#else
     double *C_work;
     double *dC_work;
#endif
    
    fp = fopen ("results_dsymv.txt", "w") ;
    if( fp == NULL ){ printf("Couldn't open output file\n"); exit(1);}

    printf("HEMV double Precision\n\n"
           "Usage\n\t\t testing_dsymv U|L N\n\n");

    N = 8*1024+64;
    if( argc > 1 ) {
      uplo = argv[1][0];
    }
    if( argc > 2 ) {
      istart = N = atoi( argv[2] );
    }
    LDA = ((N+31)/32)*32;
    matsize = N*LDA;
    vecsize = N*incx;

    TESTING_MALLOC( A, double, matsize );
    TESTING_MALLOC( X, double, vecsize );
    TESTING_MALLOC( Y, double, vecsize );
    TESTING_MALLOC( Ycublas, double, vecsize );
    TESTING_MALLOC( Ymagma,  double, vecsize );

    TESTING_DEVALLOC( dA, double, matsize );
    TESTING_DEVALLOC( dX, double, vecsize );
    TESTING_DEVALLOC( dY, double, vecsize );

#if defined(PRECISION_z) || defined(PRECISION_c)

#else
    int blocks    = N / 64 + (N % 64 != 0);
    int workspace = LDA * (blocks + 1);
    TESTING_MALLOC(    C_work, double, workspace );
    TESTING_DEVALLOC( dC_work, double, workspace );
#endif        

    /* Initialize the matrix */
    lapackf77_dlarnv( &ione, ISEED, &matsize, A );
    /* Make A hermitian */
    { 
        magma_int_t i, j;
        for(i=0; i<N; i++) {
            A[i*LDA+i] = MAGMA_D_MAKE( MAGMA_D_REAL(A[i*LDA+i]), 0. );
            for(j=0; j<i; j++)
                A[i*LDA+j] = (A[j*LDA+i]);
        }
    }
	
    printf( "   n   CUBLAS,Gflop/s   MAGMABLAS,Gflop/s      \"error\"\n" 
            "==============================================================\n");
    fprintf(fp, "   n   CUBLAS,Gflop/s   MAGMABLAS,Gflop/s      \"error\"\n" 
            "==============================================================\n");
    
    for( i = istart; i<N+1; i = (int)((i+1)*1.1) )
    {
        m = i;
	lda = ((m+31)/32)*32;
        flops = FLOPS( (double)m ) / 1e6;

        printf(      "%5d ", m );
        fprintf( fp, "%5d ", m );

        vecsize = m * incx;
        lapackf77_dlarnv( &ione, ISEED, &vecsize, X );
        lapackf77_dlarnv( &ione, ISEED, &vecsize, Y );

        /* =====================================================================
           Performs operation using CUDA-BLAS
           =================================================================== */
        cublasSetMatrix( m, m, sizeof( double ), A, LDA ,  dA, lda  );
        cublasSetVector( m,    sizeof( double ), X, incx, dX, incx );
        cublasSetVector( m,    sizeof( double ), Y, incx, dY, incx );

#if defined(PRECISION_z) || defined(PRECISION_c)

#else
    	blocks    = m / 64 + (m % 64 != 0);
    	cublasSetMatrix(lda,blocks, sizeof( double ), C_work, LDA , dC_work, lda);
#endif        
	start = get_current_time();
        cublasDsymv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        end = get_current_time();

        cublasGetVector( m, sizeof( double ), dY, incx, Ycublas, incx );
	
        
        cuda_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", cuda_perf );
        fprintf(fp, "%11.2f", cuda_perf );
	
        cublasSetVector( m, sizeof( double ), Y, incx, dY, incx );
        magmablas_dsymv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        cublasSetVector( m, sizeof( double ), Y, incx, dY, incx );
	
        
        start = get_current_time();
        magmablas_dsymv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        end = get_current_time();
        
        cublasGetVector( m, sizeof( double ), dY, incx, Ymagma, incx );

        magma_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", magma_perf );
        fprintf(fp, "%11.2f", magma_perf );

        /* =====================================================================
           Computing the Difference Cublas VS Magma
           =================================================================== */
        
        blasf77_daxpy( &m, &mzone, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_dlange( "M", &m, &ione, Ycublas, &m, work );
            
#if 0
        printf(      "\t\t %8.6e", error / m );
        fprintf( fp, "\t\t %8.6e", error / m );

        /*
         * Extra check with cblas vs magma
         */
        cblas_dcopy( m, Y, incx, Ycublas, incx );
        cblas_dsymv( CblasColMajor, CblasLower, m, 
                     (alpha), A, LDA, X, incx, 
                     (beta), Ycublas, incx );
 
        blasf77_daxpy( &m, &mzone, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_dlange( "M", &m, &ione, Ycublas, &m, work );
#endif

        printf(      "\t\t %8.6e\n", error / m );
        fprintf( fp, "\t\t %8.6e\n", error / m );
    }
    
    fclose( fp ) ; 

    /* Free Memory */
    TESTING_FREE( A );
    TESTING_FREE( X );
    TESTING_FREE( Y );
    TESTING_FREE( Ycublas );
    TESTING_FREE( Ymagma );

    TESTING_DEVFREE( dA );
    TESTING_DEVFREE( dX );
    TESTING_DEVFREE( dY );
#if defined(PRECISION_z) || defined(PRECISION_c)
	
#else 
    TESTING_FREE( C_work );
    TESTING_DEVFREE( dC_work );
#endif        

    /* Free device */
    TESTING_CUDA_FINALIZE();
    return 0;
}	
