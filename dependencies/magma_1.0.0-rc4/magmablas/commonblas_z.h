/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @precisions normal z -> s d c
 */

#ifndef _COMMONBLAS_Z_H_
#define _COMMONBLAS_Z_H_

#ifdef __cplusplus
extern "C" {
#endif

void magmablas_zgemv_MLU( magma_int_t, magma_int_t, cuDoubleComplex *, magma_int_t, cuDoubleComplex *, cuDoubleComplex *);
void magmablas_zgemv32_tesla(char, magma_int_t, magma_int_t, cuDoubleComplex, cuDoubleComplex *, magma_int_t, cuDoubleComplex *, cuDoubleComplex *);
void magmablas_zgemvt1_tesla(magma_int_t,magma_int_t,cuDoubleComplex,cuDoubleComplex *, magma_int_t,cuDoubleComplex *,cuDoubleComplex *);
void magmablas_zgemvt2_tesla(magma_int_t,magma_int_t,cuDoubleComplex,cuDoubleComplex *, magma_int_t,cuDoubleComplex *,cuDoubleComplex *);
void magmablas_zgemvt_tesla(magma_int_t,magma_int_t,cuDoubleComplex,cuDoubleComplex *, magma_int_t,cuDoubleComplex *,cuDoubleComplex *);
void magmablas_zgemv_tesla(magma_int_t M, magma_int_t N, cuDoubleComplex *A, magma_int_t lda, cuDoubleComplex *X, cuDoubleComplex *);

void magmablas_zsymv6(char, magma_int_t, cuDoubleComplex, cuDoubleComplex *, magma_int_t, cuDoubleComplex *, magma_int_t, cuDoubleComplex, cuDoubleComplex *, magma_int_t, cuDoubleComplex *, magma_int_t);

#define MAGMABLAS_ZGEMM( name ) void magmablas_zgemm_kernel_##name(cuDoubleComplex *C, const cuDoubleComplex *A, const cuDoubleComplex *B, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t lda, magma_int_t ldb, magma_int_t ldc, cuDoubleComplex alpha, cuDoubleComplex beta)
MAGMABLAS_ZGEMM( a_0                       );
MAGMABLAS_ZGEMM( ab_0                      );
MAGMABLAS_ZGEMM( N_N_64_16_16_16_4         );
MAGMABLAS_ZGEMM( N_N_64_16_16_16_4_special );
MAGMABLAS_ZGEMM( N_T_64_16_4_16_4          );
MAGMABLAS_ZGEMM( T_N_32_32_8_8_8           );
MAGMABLAS_ZGEMM( T_T_64_16_16_16_4         );
MAGMABLAS_ZGEMM( T_T_64_16_16_16_4_v2      );

void magmablas_zgemm_fermi80(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, cuDoubleComplex alpha, const cuDoubleComplex *A, magma_int_t lda, const cuDoubleComplex *B, magma_int_t ldb, cuDoubleComplex beta, cuDoubleComplex *C, magma_int_t ldc);
void magmablas_zgemm_fermi64(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, cuDoubleComplex alpha, const cuDoubleComplex *A, magma_int_t lda, const cuDoubleComplex *B, magma_int_t ldb, cuDoubleComplex beta, cuDoubleComplex *C, magma_int_t ldc);

#ifdef __cplusplus
}
#endif

#endif /* _COMMONBLAS_Z_H_ */
