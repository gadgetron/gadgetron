/*
 *  -- MAGMA (version 1.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2010
 *
 *  @generated s
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

#define PRECISION_s
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_SYMV(n) + 2. * FADDS_SYMV(n))
#else
#define FLOPS(n) (      FMULS_SYMV(n) +      FADDS_SYMV(n))
#endif

int main(int argc, char **argv)
{	
    TESTING_CUDA_INIT();

    TimeStruct  start, end;
    float      flops, magma_perf, cuda_perf, error, work[1];
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    float mzone = MAGMA_S_NEG_ONE;

    FILE        *fp ; 
    magma_int_t N, m, i, lda, LDA;
    magma_int_t matsize;
    magma_int_t vecsize;
    magma_int_t istart = 64;
    magma_int_t incx = 1;
    char        uplo = MagmaLower;
    float alpha = MAGMA_S_MAKE(1., 0.); // MAGMA_S_MAKE(  1.5, -2.3 );
    float beta  = MAGMA_S_MAKE(0., 0.); // MAGMA_S_MAKE( -0.6,  0.8 );
    float *A, *X, *Y, *Ycublas, *Ymagma;
    float *dA, *dX, *dY;
#if defined(PRECISION_z) || defined(PRECISION_c)

#else
     float *C_work;
     float *dC_work;
#endif
    
    fp = fopen ("results_ssymv.txt", "w") ;
    if( fp == NULL ){ printf("Couldn't open output file\n"); exit(1);}

    printf("HEMV float Precision\n\n"
           "Usage\n\t\t testing_ssymv U|L N\n\n");

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

    TESTING_MALLOC( A, float, matsize );
    TESTING_MALLOC( X, float, vecsize );
    TESTING_MALLOC( Y, float, vecsize );
    TESTING_MALLOC( Ycublas, float, vecsize );
    TESTING_MALLOC( Ymagma,  float, vecsize );

    TESTING_DEVALLOC( dA, float, matsize );
    TESTING_DEVALLOC( dX, float, vecsize );
    TESTING_DEVALLOC( dY, float, vecsize );

#if defined(PRECISION_z) || defined(PRECISION_c)

#else
    int blocks    = N / 64 + (N % 64 != 0);
    int workspace = LDA * (blocks + 1);
    TESTING_MALLOC(    C_work, float, workspace );
    TESTING_DEVALLOC( dC_work, float, workspace );
#endif        

    /* Initialize the matrix */
    lapackf77_slarnv( &ione, ISEED, &matsize, A );
    /* Make A hermitian */
    { 
        magma_int_t i, j;
        for(i=0; i<N; i++) {
            A[i*LDA+i] = MAGMA_S_MAKE( MAGMA_S_REAL(A[i*LDA+i]), 0. );
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
        flops = FLOPS( (float)m ) / 1e6;

        printf(      "%5d ", m );
        fprintf( fp, "%5d ", m );

        vecsize = m * incx;
        lapackf77_slarnv( &ione, ISEED, &vecsize, X );
        lapackf77_slarnv( &ione, ISEED, &vecsize, Y );

        /* =====================================================================
           Performs operation using CUDA-BLAS
           =================================================================== */
        cublasSetMatrix( m, m, sizeof( float ), A, LDA ,  dA, lda  );
        cublasSetVector( m,    sizeof( float ), X, incx, dX, incx );
        cublasSetVector( m,    sizeof( float ), Y, incx, dY, incx );

#if defined(PRECISION_z) || defined(PRECISION_c)

#else
    	blocks    = m / 64 + (m % 64 != 0);
    	cublasSetMatrix(lda,blocks, sizeof( float ), C_work, LDA , dC_work, lda);
#endif        
	start = get_current_time();
        cublasSsymv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        end = get_current_time();

        cublasGetVector( m, sizeof( float ), dY, incx, Ycublas, incx );
	
        
        cuda_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", cuda_perf );
        fprintf(fp, "%11.2f", cuda_perf );
	
        cublasSetVector( m, sizeof( float ), Y, incx, dY, incx );
        magmablas_ssymv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        cublasSetVector( m, sizeof( float ), Y, incx, dY, incx );
	
        
        start = get_current_time();
        magmablas_ssymv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        end = get_current_time();
        
        cublasGetVector( m, sizeof( float ), dY, incx, Ymagma, incx );

        magma_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", magma_perf );
        fprintf(fp, "%11.2f", magma_perf );

        /* =====================================================================
           Computing the Difference Cublas VS Magma
           =================================================================== */
        
        blasf77_saxpy( &m, &mzone, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_slange( "M", &m, &ione, Ycublas, &m, work );
            
#if 0
        printf(      "\t\t %8.6e", error / m );
        fprintf( fp, "\t\t %8.6e", error / m );

        /*
         * Extra check with cblas vs magma
         */
        cblas_scopy( m, Y, incx, Ycublas, incx );
        cblas_ssymv( CblasColMajor, CblasLower, m, 
                     (alpha), A, LDA, X, incx, 
                     (beta), Ycublas, incx );
 
        blasf77_saxpy( &m, &mzone, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_slange( "M", &m, &ione, Ycublas, &m, work );
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
