/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated c
 */

#ifndef _COMMONBLAS_C_H_
#define _COMMONBLAS_C_H_

#ifdef __cplusplus
extern "C" {
#endif

void magmablas_cgemv_MLU( magma_int_t, magma_int_t, cuFloatComplex *, magma_int_t, cuFloatComplex *, cuFloatComplex *);
void magmablas_cgemv32_tesla(char, magma_int_t, magma_int_t, cuFloatComplex, cuFloatComplex *, magma_int_t, cuFloatComplex *, cuFloatComplex *);
void magmablas_cgemvt1_tesla(magma_int_t,magma_int_t,cuFloatComplex,cuFloatComplex *, magma_int_t,cuFloatComplex *,cuFloatComplex *);
void magmablas_cgemvt2_tesla(magma_int_t,magma_int_t,cuFloatComplex,cuFloatComplex *, magma_int_t,cuFloatComplex *,cuFloatComplex *);
void magmablas_cgemvt_tesla(magma_int_t,magma_int_t,cuFloatComplex,cuFloatComplex *, magma_int_t,cuFloatComplex *,cuFloatComplex *);
void magmablas_cgemv_tesla(magma_int_t M, magma_int_t N, cuFloatComplex *A, magma_int_t lda, cuFloatComplex *X, cuFloatComplex *);

void magmablas_csymv6(char, magma_int_t, cuFloatComplex, cuFloatComplex *, magma_int_t, cuFloatComplex *, magma_int_t, cuFloatComplex, cuFloatComplex *, magma_int_t, cuFloatComplex *, magma_int_t);

#define MAGMABLAS_CGEMM( name ) void magmablas_cgemm_kernel_##name(cuFloatComplex *C, const cuFloatComplex *A, const cuFloatComplex *B, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t lda, magma_int_t ldb, magma_int_t ldc, cuFloatComplex alpha, cuFloatComplex beta)
MAGMABLAS_CGEMM( a_0                       );
MAGMABLAS_CGEMM( ab_0                      );
MAGMABLAS_CGEMM( N_N_64_16_16_16_4         );
MAGMABLAS_CGEMM( N_N_64_16_16_16_4_special );
MAGMABLAS_CGEMM( N_T_64_16_4_16_4          );
MAGMABLAS_CGEMM( T_N_32_32_8_8_8           );
MAGMABLAS_CGEMM( T_T_64_16_16_16_4         );
MAGMABLAS_CGEMM( T_T_64_16_16_16_4_v2      );

void magmablas_cgemm_fermi80(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, cuFloatComplex alpha, const cuFloatComplex *A, magma_int_t lda, const cuFloatComplex *B, magma_int_t ldb, cuFloatComplex beta, cuFloatComplex *C, magma_int_t ldc);
void magmablas_cgemm_fermi64(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, cuFloatComplex alpha, const cuFloatComplex *A, magma_int_t lda, const cuFloatComplex *B, magma_int_t ldb, cuFloatComplex beta, cuFloatComplex *C, magma_int_t ldc);

#ifdef __cplusplus
}
#endif

#endif /* _COMMONBLAS_C_H_ */
