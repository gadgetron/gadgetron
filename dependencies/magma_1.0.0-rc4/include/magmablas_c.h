/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated c
 */

#ifndef _MAGMABLAS_C_H_
#define _MAGMABLAS_C_H_

#define PRECISION_c
#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
float cpu_gpu_cdiff(             int M, int N, 
				  cuFloatComplex * a, int lda, 
				  cuFloatComplex *da, int ldda);
void czero_32x32_block(           cuFloatComplex *, magma_int_t);
void czero_nbxnb_block(           magma_int_t, cuFloatComplex *, magma_int_t);
void magmablas_cinplace_transpose(cuFloatComplex *, magma_int_t, magma_int_t);
void magmablas_cpermute_long(     cuFloatComplex *, magma_int_t, 
				  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_cpermute_long2(    cuFloatComplex *, magma_int_t, 
				  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_ctranspose(        cuFloatComplex *, magma_int_t, 
				  cuFloatComplex *, magma_int_t, 
				  magma_int_t, magma_int_t);
void magmablas_ctranspose2(       cuFloatComplex *, magma_int_t, 
				  cuFloatComplex *, magma_int_t, 
				  magma_int_t, magma_int_t);
  
  /*
   * LAPACK auxiliary functions
   */
void   magmablas_clacpy( char uplo, 
			 magma_int_t m, magma_int_t n, 
			 cuFloatComplex *A, magma_int_t lda, 
			 cuFloatComplex *B, magma_int_t ldb);
float magmablas_clange( char norm, 
			 magma_int_t m, magma_int_t n, 
			 cuFloatComplex *A, magma_int_t lda, float *WORK);
float magmablas_clanhe( char norm, char uplo, 
			 magma_int_t n,
			 cuFloatComplex *A, magma_int_t lda, float *WORK);
float magmablas_clansy( char norm, char uplo,
			 magma_int_t n, 
			 cuFloatComplex *A, magma_int_t lda, float *WORK);
void   magmablas_claset( magma_int_t m, magma_int_t n,
			 cuFloatComplex *A, magma_int_t lda);
void   magmablas_claswp( magma_int_t N, 
			 cuFloatComplex *dAT, magma_int_t lda, 
			 magma_int_t i1,  magma_int_t i2, 
			 magma_int_t *ipiv, magma_int_t inci );
void   magmablas_claswpx(magma_int_t N, 
			 cuFloatComplex *dAT, magma_int_t ldx, magma_int_t ldy, 
			 magma_int_t i1, magma_int_t i2,
			 magma_int_t *ipiv, magma_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magmablas_cswap(   magma_int_t N, 
			  cuFloatComplex *dA1, magma_int_t lda1, 
			  cuFloatComplex *dA2, magma_int_t lda2 );
void   magmablas_cswapblk(char storev, 
			  magma_int_t N, 
			  cuFloatComplex *dA1, magma_int_t lda1, 
			  cuFloatComplex *dA2, magma_int_t lda2,
			  magma_int_t i1, magma_int_t i2, 
			  magma_int_t *ipiv, magma_int_t inci, 
			  magma_int_t offset);

  /*
   * Level 2 BLAS
   */
void magmablas_cgemv(char t, magma_int_t M, magma_int_t N, 
		     cuFloatComplex alpha,
		     cuFloatComplex *A, magma_int_t lda, 
		     cuFloatComplex * X, magma_int_t incX, 
		     cuFloatComplex beta, 
		     cuFloatComplex *Y, magma_int_t incY);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magmablas_chemv(char u, magma_int_t N, 
                            cuFloatComplex alpha, 
                            cuFloatComplex *A, magma_int_t lda, 
                            cuFloatComplex *X, magma_int_t incX, 
                            cuFloatComplex beta, 
                            cuFloatComplex *Y, magma_int_t incY);
#endif
magma_int_t magmablas_csymv(char u, magma_int_t N, 
                            cuFloatComplex alpha, 
                            cuFloatComplex *A, magma_int_t lda, 
                            cuFloatComplex *X, magma_int_t incX, 
                            cuFloatComplex beta, 
                            cuFloatComplex *Y, magma_int_t incY);

  /*
   * Level 3 BLAS
   */
void magmablas_cgemm(char tA, char tB,
		     magma_int_t m, magma_int_t n, magma_int_t k, 
		     cuFloatComplex alpha,
		     const cuFloatComplex *A, magma_int_t lda, 
		     const cuFloatComplex *B, magma_int_t ldb, 
		     cuFloatComplex beta,
		     cuFloatComplex *C, magma_int_t ldc);
void magmablas_cgemm_fermi80(char tA, char tB, 
			     magma_int_t m, magma_int_t n, magma_int_t k,
			     cuFloatComplex alpha, 
			     const cuFloatComplex *A, magma_int_t lda, 
			     const cuFloatComplex *B, magma_int_t ldb,
			     cuFloatComplex beta, 
			     cuFloatComplex *C, magma_int_t ldc);
void magmablas_cgemm_fermi64(char tA, char tB, 
			     magma_int_t m, magma_int_t n, magma_int_t k,
			     cuFloatComplex alpha, 
			     const cuFloatComplex *A, magma_int_t lda, 
			     const cuFloatComplex *B, magma_int_t ldb, 
			     cuFloatComplex beta, 
			     cuFloatComplex *C, magma_int_t ldc);
void magmablas_chemm(char s, char u,          
		     magma_int_t m, magma_int_t n,
		     cuFloatComplex alpha, 
		     const cuFloatComplex *A, magma_int_t lda,
		     const cuFloatComplex *B, magma_int_t ldb,
		     cuFloatComplex beta, 
		     cuFloatComplex *C, magma_int_t ldc);
void magmablas_csymm(char s, char u,
		     magma_int_t m, magma_int_t n,
		     cuFloatComplex alpha, 
		     const cuFloatComplex *A, magma_int_t lda, 
		     const cuFloatComplex *B, magma_int_t ldb,
		     cuFloatComplex beta,
		     cuFloatComplex *C, magma_int_t ldc);
void magmablas_csyrk(char u, char t,
		     magma_int_t n, magma_int_t k, 
		     cuFloatComplex alpha, 
		     const cuFloatComplex *A, magma_int_t lda,
		     cuFloatComplex beta,
		     cuFloatComplex *C, magma_int_t ldc);
void magmablas_cherk(char u, char t,
		     magma_int_t n, magma_int_t k, 
		     float  alpha, 
		     const cuFloatComplex *A, magma_int_t lda,
		     float  beta, 
		     cuFloatComplex *C, magma_int_t ldc);
void magmablas_csyr2k(char u, char t,
		      magma_int_t n, magma_int_t k,
		      cuFloatComplex alpha, 
		      const cuFloatComplex *A, magma_int_t lda,
		      const cuFloatComplex *B, magma_int_t ldb, 
		      cuFloatComplex beta, 
		      cuFloatComplex *C, magma_int_t ldc);
void magmablas_cher2k(char u, char t,
		      magma_int_t n, magma_int_t k, 
		      cuFloatComplex alpha, 
		      const cuFloatComplex *A, magma_int_t lda, 
		      const cuFloatComplex *B, magma_int_t ldb,
		      float  beta,
		      cuFloatComplex *C, magma_int_t ldc);
void magmablas_ctrmm(char s, char u, char t,  char d, 
		     magma_int_t m, magma_int_t n,
		     cuFloatComplex alpha,
		     const cuFloatComplex *A, magma_int_t lda,
		     cuFloatComplex *B, magma_int_t ldb);
void magmablas_ctrsm(char s, char u, char t, char d,
		     magma_int_t m, magma_int_t n,
		     cuFloatComplex alpha,
		     /*const*/ cuFloatComplex *A, magma_int_t lda,
		     cuFloatComplex *B, magma_int_t ldb);


  /*
   * Workspace interface (alphabetical order)
   */
magma_int_t magmablasw_csymv(char u, magma_int_t N, 
			     cuFloatComplex alpha, 
			     cuFloatComplex *A, magma_int_t lda, 
			     cuFloatComplex *X, magma_int_t incX, 
			     cuFloatComplex beta, 
			     cuFloatComplex *Y, magma_int_t incY,
			     cuFloatComplex *dWork);

#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif
