/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated s
 */

#ifndef _COMMONBLAS_S_H_
#define _COMMONBLAS_S_H_

#ifdef __cplusplus
extern "C" {
#endif

void magmablas_sgemv_MLU( magma_int_t, magma_int_t, float *, magma_int_t, float *, float *);
void magmablas_sgemv32_tesla(char, magma_int_t, magma_int_t, float, float *, magma_int_t, float *, float *);
void magmablas_sgemvt1_tesla(magma_int_t,magma_int_t,float,float *, magma_int_t,float *,float *);
void magmablas_sgemvt2_tesla(magma_int_t,magma_int_t,float,float *, magma_int_t,float *,float *);
void magmablas_sgemvt_tesla(magma_int_t,magma_int_t,float,float *, magma_int_t,float *,float *);
void magmablas_sgemv_tesla(magma_int_t M, magma_int_t N, float *A, magma_int_t lda, float *X, float *);

void magmablas_ssymv6(char, magma_int_t, float, float *, magma_int_t, float *, magma_int_t, float, float *, magma_int_t, float *, magma_int_t);

#define MAGMABLAS_SGEMM( name ) void magmablas_sgemm_kernel_##name(float *C, const float *A, const float *B, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t lda, magma_int_t ldb, magma_int_t ldc, float alpha, float beta)
MAGMABLAS_SGEMM( a_0                       );
MAGMABLAS_SGEMM( ab_0                      );
MAGMABLAS_SGEMM( N_N_64_16_16_16_4         );
MAGMABLAS_SGEMM( N_N_64_16_16_16_4_special );
MAGMABLAS_SGEMM( N_T_64_16_4_16_4          );
MAGMABLAS_SGEMM( T_N_32_32_8_8_8           );
MAGMABLAS_SGEMM( T_T_64_16_16_16_4         );
MAGMABLAS_SGEMM( T_T_64_16_16_16_4_v2      );

void magmablas_sgemm_fermi80(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, float alpha, const float *A, magma_int_t lda, const float *B, magma_int_t ldb, float beta, float *C, magma_int_t ldc);
void magmablas_sgemm_fermi64(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, float alpha, const float *A, magma_int_t lda, const float *B, magma_int_t ldb, float beta, float *C, magma_int_t ldc);

#ifdef __cplusplus
}
#endif

#endif /* _COMMONBLAS_S_H_ */
