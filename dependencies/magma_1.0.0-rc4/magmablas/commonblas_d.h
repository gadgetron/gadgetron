/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated d
 */

#ifndef _COMMONBLAS_D_H_
#define _COMMONBLAS_D_H_

#ifdef __cplusplus
extern "C" {
#endif

void magmablas_dgemv_MLU( magma_int_t, magma_int_t, double *, magma_int_t, double *, double *);
void magmablas_dgemv32_tesla(char, magma_int_t, magma_int_t, double, double *, magma_int_t, double *, double *);
void magmablas_dgemvt1_tesla(magma_int_t,magma_int_t,double,double *, magma_int_t,double *,double *);
void magmablas_dgemvt2_tesla(magma_int_t,magma_int_t,double,double *, magma_int_t,double *,double *);
void magmablas_dgemvt_tesla(magma_int_t,magma_int_t,double,double *, magma_int_t,double *,double *);
void magmablas_dgemv_tesla(magma_int_t M, magma_int_t N, double *A, magma_int_t lda, double *X, double *);

void magmablas_dsymv6(char, magma_int_t, double, double *, magma_int_t, double *, magma_int_t, double, double *, magma_int_t, double *, magma_int_t);

#define MAGMABLAS_DGEMM( name ) void magmablas_dgemm_kernel_##name(double *C, const double *A, const double *B, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t lda, magma_int_t ldb, magma_int_t ldc, double alpha, double beta)
MAGMABLAS_DGEMM( a_0                       );
MAGMABLAS_DGEMM( ab_0                      );
MAGMABLAS_DGEMM( N_N_64_16_16_16_4         );
MAGMABLAS_DGEMM( N_N_64_16_16_16_4_special );
MAGMABLAS_DGEMM( N_T_64_16_4_16_4          );
MAGMABLAS_DGEMM( T_N_32_32_8_8_8           );
MAGMABLAS_DGEMM( T_T_64_16_16_16_4         );
MAGMABLAS_DGEMM( T_T_64_16_16_16_4_v2      );

void magmablas_dgemm_fermi80(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, double alpha, const double *A, magma_int_t lda, const double *B, magma_int_t ldb, double beta, double *C, magma_int_t ldc);
void magmablas_dgemm_fermi64(char tA, char tB, magma_int_t m, magma_int_t n, magma_int_t k, double alpha, const double *A, magma_int_t lda, const double *B, magma_int_t ldb, double beta, double *C, magma_int_t ldc);

#ifdef __cplusplus
}
#endif

#endif /* _COMMONBLAS_D_H_ */
