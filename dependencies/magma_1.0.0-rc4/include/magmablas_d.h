/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated d
 */

#ifndef _MAGMABLAS_D_H_
#define _MAGMABLAS_D_H_

#define PRECISION_d
#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
double cpu_gpu_ddiff(             int M, int N, 
				  double * a, int lda, 
				  double *da, int ldda);
void dzero_32x32_block(           double *, magma_int_t);
void dzero_nbxnb_block(           magma_int_t, double *, magma_int_t);
void magmablas_dinplace_transpose(double *, magma_int_t, magma_int_t);
void magmablas_dpermute_long(     double *, magma_int_t, 
				  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_dpermute_long2(    double *, magma_int_t, 
				  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_dtranspose(        double *, magma_int_t, 
				  double *, magma_int_t, 
				  magma_int_t, magma_int_t);
void magmablas_dtranspose2(       double *, magma_int_t, 
				  double *, magma_int_t, 
				  magma_int_t, magma_int_t);
  
  /*
   * LAPACK auxiliary functions
   */
void   magmablas_dlacpy( char uplo, 
			 magma_int_t m, magma_int_t n, 
			 double *A, magma_int_t lda, 
			 double *B, magma_int_t ldb);
double magmablas_dlange( char norm, 
			 magma_int_t m, magma_int_t n, 
			 double *A, magma_int_t lda, double *WORK);
double magmablas_dlansy( char norm, char uplo, 
			 magma_int_t n,
			 double *A, magma_int_t lda, double *WORK);
double magmablas_dlansy( char norm, char uplo,
			 magma_int_t n, 
			 double *A, magma_int_t lda, double *WORK);
void   magmablas_dlaset( magma_int_t m, magma_int_t n,
			 double *A, magma_int_t lda);
void   magmablas_dlaswp( magma_int_t N, 
			 double *dAT, magma_int_t lda, 
			 magma_int_t i1,  magma_int_t i2, 
			 magma_int_t *ipiv, magma_int_t inci );
void   magmablas_dlaswpx(magma_int_t N, 
			 double *dAT, magma_int_t ldx, magma_int_t ldy, 
			 magma_int_t i1, magma_int_t i2,
			 magma_int_t *ipiv, magma_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magmablas_dswap(   magma_int_t N, 
			  double *dA1, magma_int_t lda1, 
			  double *dA2, magma_int_t lda2 );
void   magmablas_dswapblk(char storev, 
			  magma_int_t N, 
			  double *dA1, magma_int_t lda1, 
			  double *dA2, magma_int_t lda2,
			  magma_int_t i1, magma_int_t i2, 
			  magma_int_t *ipiv, magma_int_t inci, 
			  magma_int_t offset);

  /*
   * Level 2 BLAS
   */
void magmablas_dgemv(char t, magma_int_t M, magma_int_t N, 
		     double alpha,
		     double *A, magma_int_t lda, 
		     double * X, magma_int_t incX, 
		     double beta, 
		     double *Y, magma_int_t incY);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magmablas_dsymv(char u, magma_int_t N, 
                            double alpha, 
                            double *A, magma_int_t lda, 
                            double *X, magma_int_t incX, 
                            double beta, 
                            double *Y, magma_int_t incY);
#endif
magma_int_t magmablas_dsymv(char u, magma_int_t N, 
                            double alpha, 
                            double *A, magma_int_t lda, 
                            double *X, magma_int_t incX, 
                            double beta, 
                            double *Y, magma_int_t incY);

  /*
   * Level 3 BLAS
   */
void magmablas_dgemm(char tA, char tB,
		     magma_int_t m, magma_int_t n, magma_int_t k, 
		     double alpha,
		     const double *A, magma_int_t lda, 
		     const double *B, magma_int_t ldb, 
		     double beta,
		     double *C, magma_int_t ldc);
void magmablas_dgemm_fermi80(char tA, char tB, 
			     magma_int_t m, magma_int_t n, magma_int_t k,
			     double alpha, 
			     const double *A, magma_int_t lda, 
			     const double *B, magma_int_t ldb,
			     double beta, 
			     double *C, magma_int_t ldc);
void magmablas_dgemm_fermi64(char tA, char tB, 
			     magma_int_t m, magma_int_t n, magma_int_t k,
			     double alpha, 
			     const double *A, magma_int_t lda, 
			     const double *B, magma_int_t ldb, 
			     double beta, 
			     double *C, magma_int_t ldc);
void magmablas_dsymm(char s, char u,          
		     magma_int_t m, magma_int_t n,
		     double alpha, 
		     const double *A, magma_int_t lda,
		     const double *B, magma_int_t ldb,
		     double beta, 
		     double *C, magma_int_t ldc);
void magmablas_dsymm(char s, char u,
		     magma_int_t m, magma_int_t n,
		     double alpha, 
		     const double *A, magma_int_t lda, 
		     const double *B, magma_int_t ldb,
		     double beta,
		     double *C, magma_int_t ldc);
void magmablas_dsyrk(char u, char t,
		     magma_int_t n, magma_int_t k, 
		     double alpha, 
		     const double *A, magma_int_t lda,
		     double beta,
		     double *C, magma_int_t ldc);
void magmablas_dsyrk(char u, char t,
		     magma_int_t n, magma_int_t k, 
		     double  alpha, 
		     const double *A, magma_int_t lda,
		     double  beta, 
		     double *C, magma_int_t ldc);
void magmablas_dsyr2k(char u, char t,
		      magma_int_t n, magma_int_t k,
		      double alpha, 
		      const double *A, magma_int_t lda,
		      const double *B, magma_int_t ldb, 
		      double beta, 
		      double *C, magma_int_t ldc);
void magmablas_dsyr2k(char u, char t,
		      magma_int_t n, magma_int_t k, 
		      double alpha, 
		      const double *A, magma_int_t lda, 
		      const double *B, magma_int_t ldb,
		      double  beta,
		      double *C, magma_int_t ldc);
void magmablas_dtrmm(char s, char u, char t,  char d, 
		     magma_int_t m, magma_int_t n,
		     double alpha,
		     const double *A, magma_int_t lda,
		     double *B, magma_int_t ldb);
void magmablas_dtrsm(char s, char u, char t, char d,
		     magma_int_t m, magma_int_t n,
		     double alpha,
		     /*const*/ double *A, magma_int_t lda,
		     double *B, magma_int_t ldb);


  /*
   * Workspace interface (alphabetical order)
   */
magma_int_t magmablasw_dsymv(char u, magma_int_t N, 
			     double alpha, 
			     double *A, magma_int_t lda, 
			     double *X, magma_int_t incX, 
			     double beta, 
			     double *Y, magma_int_t incY,
			     double *dWork);

#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif
