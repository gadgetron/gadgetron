/*
 *  -- MAGMA (version 1.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2010
 *
 *  @precisions normal z -> c d s
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

#define PRECISION_z
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
    cuDoubleComplex mzone = MAGMA_Z_NEG_ONE;

    FILE        *fp ; 
    magma_int_t N, m, i, lda, LDA;
    magma_int_t matsize;
    magma_int_t vecsize;
    magma_int_t istart = 64;
    magma_int_t incx = 1;
    char        uplo = MagmaLower;
    cuDoubleComplex alpha = MAGMA_Z_MAKE(1., 0.); // MAGMA_Z_MAKE(  1.5, -2.3 );
    cuDoubleComplex beta  = MAGMA_Z_MAKE(0., 0.); // MAGMA_Z_MAKE( -0.6,  0.8 );
    cuDoubleComplex *A, *X, *Y, *Ycublas, *Ymagma;
    cuDoubleComplex *dA, *dX, *dY;
#if defined(PRECISION_z) || defined(PRECISION_c)

#else
     cuDoubleComplex *C_work;
     cuDoubleComplex *dC_work;
#endif
    
    fp = fopen ("results_zhemv.txt", "w") ;
    if( fp == NULL ){ printf("Couldn't open output file\n"); exit(1);}

    printf("HEMV cuDoubleComplex Precision\n\n"
           "Usage\n\t\t testing_zhemv U|L N\n\n");

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

    TESTING_MALLOC( A, cuDoubleComplex, matsize );
    TESTING_MALLOC( X, cuDoubleComplex, vecsize );
    TESTING_MALLOC( Y, cuDoubleComplex, vecsize );
    TESTING_MALLOC( Ycublas, cuDoubleComplex, vecsize );
    TESTING_MALLOC( Ymagma,  cuDoubleComplex, vecsize );

    TESTING_DEVALLOC( dA, cuDoubleComplex, matsize );
    TESTING_DEVALLOC( dX, cuDoubleComplex, vecsize );
    TESTING_DEVALLOC( dY, cuDoubleComplex, vecsize );

#if defined(PRECISION_z) || defined(PRECISION_c)

#else
    int blocks    = N / 64 + (N % 64 != 0);
    int workspace = LDA * (blocks + 1);
    TESTING_MALLOC(    C_work, cuDoubleComplex, workspace );
    TESTING_DEVALLOC( dC_work, cuDoubleComplex, workspace );
#endif        

    /* Initialize the matrix */
    lapackf77_zlarnv( &ione, ISEED, &matsize, A );
    /* Make A hermitian */
    { 
        magma_int_t i, j;
        for(i=0; i<N; i++) {
            A[i*LDA+i] = MAGMA_Z_MAKE( MAGMA_Z_REAL(A[i*LDA+i]), 0. );
            for(j=0; j<i; j++)
                A[i*LDA+j] = cuConj(A[j*LDA+i]);
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
        lapackf77_zlarnv( &ione, ISEED, &vecsize, X );
        lapackf77_zlarnv( &ione, ISEED, &vecsize, Y );

        /* =====================================================================
           Performs operation using CUDA-BLAS
           =================================================================== */
        cublasSetMatrix( m, m, sizeof( cuDoubleComplex ), A, LDA ,  dA, lda  );
        cublasSetVector( m,    sizeof( cuDoubleComplex ), X, incx, dX, incx );
        cublasSetVector( m,    sizeof( cuDoubleComplex ), Y, incx, dY, incx );

#if defined(PRECISION_z) || defined(PRECISION_c)

#else
    	blocks    = m / 64 + (m % 64 != 0);
    	cublasSetMatrix(lda,blocks, sizeof( cuDoubleComplex ), C_work, LDA , dC_work, lda);
#endif        
	start = get_current_time();
        cublasZhemv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        end = get_current_time();

        cublasGetVector( m, sizeof( cuDoubleComplex ), dY, incx, Ycublas, incx );
	
        
        cuda_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", cuda_perf );
        fprintf(fp, "%11.2f", cuda_perf );
	
        cublasSetVector( m, sizeof( cuDoubleComplex ), Y, incx, dY, incx );
        magmablas_zhemv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        cublasSetVector( m, sizeof( cuDoubleComplex ), Y, incx, dY, incx );
	
        
        start = get_current_time();
        magmablas_zhemv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        end = get_current_time();
        
        cublasGetVector( m, sizeof( cuDoubleComplex ), dY, incx, Ymagma, incx );

        magma_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", magma_perf );
        fprintf(fp, "%11.2f", magma_perf );

        /* =====================================================================
           Computing the Difference Cublas VS Magma
           =================================================================== */
        
        blasf77_zaxpy( &m, &mzone, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_zlange( "M", &m, &ione, Ycublas, &m, work );
            
#if 0
        printf(      "\t\t %8.6e", error / m );
        fprintf( fp, "\t\t %8.6e", error / m );

        /*
         * Extra check with cblas vs magma
         */
        cblas_zcopy( m, Y, incx, Ycublas, incx );
        cblas_zhemv( CblasColMajor, CblasLower, m, 
                     CBLAS_SADDR(alpha), A, LDA, X, incx, 
                     CBLAS_SADDR(beta), Ycublas, incx );
 
        blasf77_zaxpy( &m, &mzone, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_zlange( "M", &m, &ione, Ycublas, &m, work );
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
