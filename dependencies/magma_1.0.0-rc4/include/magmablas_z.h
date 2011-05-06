/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMABLAS_Z_H_
#define _MAGMABLAS_Z_H_

#define PRECISION_z
#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
double cpu_gpu_zdiff(             int M, int N, 
				  cuDoubleComplex * a, int lda, 
				  cuDoubleComplex *da, int ldda);
void zzero_32x32_block(           cuDoubleComplex *, magma_int_t);
void zzero_nbxnb_block(           magma_int_t, cuDoubleComplex *, magma_int_t);
void magmablas_zinplace_transpose(cuDoubleComplex *, magma_int_t, magma_int_t);
void magmablas_zpermute_long(     cuDoubleComplex *, magma_int_t, 
				  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_zpermute_long2(    cuDoubleComplex *, magma_int_t, 
				  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_ztranspose(        cuDoubleComplex *, magma_int_t, 
				  cuDoubleComplex *, magma_int_t, 
				  magma_int_t, magma_int_t);
void magmablas_ztranspose2(       cuDoubleComplex *, magma_int_t, 
				  cuDoubleComplex *, magma_int_t, 
				  magma_int_t, magma_int_t);
  
  /*
   * LAPACK auxiliary functions
   */
void   magmablas_zlacpy( char uplo, 
			 magma_int_t m, magma_int_t n, 
			 cuDoubleComplex *A, magma_int_t lda, 
			 cuDoubleComplex *B, magma_int_t ldb);
double magmablas_zlange( char norm, 
			 magma_int_t m, magma_int_t n, 
			 cuDoubleComplex *A, magma_int_t lda, double *WORK);
double magmablas_zlanhe( char norm, char uplo, 
			 magma_int_t n,
			 cuDoubleComplex *A, magma_int_t lda, double *WORK);
double magmablas_zlansy( char norm, char uplo,
			 magma_int_t n, 
			 cuDoubleComplex *A, magma_int_t lda, double *WORK);
void   magmablas_zlaset( magma_int_t m, magma_int_t n,
			 cuDoubleComplex *A, magma_int_t lda);
void   magmablas_zlaswp( magma_int_t N, 
			 cuDoubleComplex *dAT, magma_int_t lda, 
			 magma_int_t i1,  magma_int_t i2, 
			 magma_int_t *ipiv, magma_int_t inci );
void   magmablas_zlaswpx(magma_int_t N, 
			 cuDoubleComplex *dAT, magma_int_t ldx, magma_int_t ldy, 
			 magma_int_t i1, magma_int_t i2,
			 magma_int_t *ipiv, magma_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magmablas_zswap(   magma_int_t N, 
			  cuDoubleComplex *dA1, magma_int_t lda1, 
			  cuDoubleComplex *dA2, magma_int_t lda2 );
void   magmablas_zswapblk(char storev, 
			  magma_int_t N, 
			  cuDoubleComplex *dA1, magma_int_t lda1, 
			  cuDoubleComplex *dA2, magma_int_t lda2,
			  magma_int_t i1, magma_int_t i2, 
			  magma_int_t *ipiv, magma_int_t inci, 
			  magma_int_t offset);

  /*
   * Level 2 BLAS
   */
void magmablas_zgemv(char t, magma_int_t M, magma_int_t N, 
		     cuDoubleComplex alpha,
		     cuDoubleComplex *A, magma_int_t lda, 
		     cuDoubleComplex * X, magma_int_t incX, 
		     cuDoubleComplex beta, 
		     cuDoubleComplex *Y, magma_int_t incY);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magmablas_zhemv(char u, magma_int_t N, 
                            cuDoubleComplex alpha, 
                            cuDoubleComplex *A, magma_int_t lda, 
                            cuDoubleComplex *X, magma_int_t incX, 
                            cuDoubleComplex beta, 
                            cuDoubleComplex *Y, magma_int_t incY);
#endif
magma_int_t magmablas_zsymv(char u, magma_int_t N, 
                            cuDoubleComplex alpha, 
                            cuDoubleComplex *A, magma_int_t lda, 
                            cuDoubleComplex *X, magma_int_t incX, 
                            cuDoubleComplex beta, 
                            cuDoubleComplex *Y, magma_int_t incY);

  /*
   * Level 3 BLAS
   */
void magmablas_zgemm(char tA, char tB,
		     magma_int_t m, magma_int_t n, magma_int_t k, 
		     cuDoubleComplex alpha,
		     const cuDoubleComplex *A, magma_int_t lda, 
		     const cuDoubleComplex *B, magma_int_t ldb, 
		     cuDoubleComplex beta,
		     cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zgemm_fermi80(char tA, char tB, 
			     magma_int_t m, magma_int_t n, magma_int_t k,
			     cuDoubleComplex alpha, 
			     const cuDoubleComplex *A, magma_int_t lda, 
			     const cuDoubleComplex *B, magma_int_t ldb,
			     cuDoubleComplex beta, 
			     cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zgemm_fermi64(char tA, char tB, 
			     magma_int_t m, magma_int_t n, magma_int_t k,
			     cuDoubleComplex alpha, 
			     const cuDoubleComplex *A, magma_int_t lda, 
			     const cuDoubleComplex *B, magma_int_t ldb, 
			     cuDoubleComplex beta, 
			     cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zhemm(char s, char u,          
		     magma_int_t m, magma_int_t n,
		     cuDoubleComplex alpha, 
		     const cuDoubleComplex *A, magma_int_t lda,
		     const cuDoubleComplex *B, magma_int_t ldb,
		     cuDoubleComplex beta, 
		     cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zsymm(char s, char u,
		     magma_int_t m, magma_int_t n,
		     cuDoubleComplex alpha, 
		     const cuDoubleComplex *A, magma_int_t lda, 
		     const cuDoubleComplex *B, magma_int_t ldb,
		     cuDoubleComplex beta,
		     cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zsyrk(char u, char t,
		     magma_int_t n, magma_int_t k, 
		     cuDoubleComplex alpha, 
		     const cuDoubleComplex *A, magma_int_t lda,
		     cuDoubleComplex beta,
		     cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zherk(char u, char t,
		     magma_int_t n, magma_int_t k, 
		     double  alpha, 
		     const cuDoubleComplex *A, magma_int_t lda,
		     double  beta, 
		     cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zsyr2k(char u, char t,
		      magma_int_t n, magma_int_t k,
		      cuDoubleComplex alpha, 
		      const cuDoubleComplex *A, magma_int_t lda,
		      const cuDoubleComplex *B, magma_int_t ldb, 
		      cuDoubleComplex beta, 
		      cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zher2k(char u, char t,
		      magma_int_t n, magma_int_t k, 
		      cuDoubleComplex alpha, 
		      const cuDoubleComplex *A, magma_int_t lda, 
		      const cuDoubleComplex *B, magma_int_t ldb,
		      double  beta,
		      cuDoubleComplex *C, magma_int_t ldc);
void magmablas_ztrmm(char s, char u, char t,  char d, 
		     magma_int_t m, magma_int_t n,
		     cuDoubleComplex alpha,
		     const cuDoubleComplex *A, magma_int_t lda,
		     cuDoubleComplex *B, magma_int_t ldb);
void magmablas_ztrsm(char s, char u, char t, char d,
		     magma_int_t m, magma_int_t n,
		     cuDoubleComplex alpha,
		     /*const*/ cuDoubleComplex *A, magma_int_t lda,
		     cuDoubleComplex *B, magma_int_t ldb);


  /*
   * Workspace interface (alphabetical order)
   */
magma_int_t magmablasw_zsymv(char u, magma_int_t N, 
			     cuDoubleComplex alpha, 
			     cuDoubleComplex *A, magma_int_t lda, 
			     cuDoubleComplex *X, magma_int_t incX, 
			     cuDoubleComplex beta, 
			     cuDoubleComplex *Y, magma_int_t incY,
			     cuDoubleComplex *dWork);

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif
