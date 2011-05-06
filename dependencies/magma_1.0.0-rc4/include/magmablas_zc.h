/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @precisions mixed zc -> ds
 */

#ifndef _MAGMABLAS_ZC_H_
#define _MAGMABLAS_ZC_H_

#include "cublas.h"
#include "cuda.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magmablas_zcaxpycp(cuFloatComplex *, cuDoubleComplex *, magma_int_t, cuDoubleComplex *, cuDoubleComplex *);
void magmablas_zaxpycp(cuDoubleComplex *, cuDoubleComplex *, magma_int_t, cuDoubleComplex *);
void magmablas_zclaswp(magma_int_t, cuDoubleComplex *, magma_int_t, cuFloatComplex *, magma_int_t, magma_int_t *, magma_int_t);
void magmablas_zlag2c(magma_int_t M, magma_int_t N, const cuDoubleComplex *A, magma_int_t lda,  cuFloatComplex *SA, magma_int_t ldsa, magma_int_t *info);

void magmablas_clag2z(magma_int_t M, magma_int_t N, 
                      cuFloatComplex  *SA, magma_int_t ldsa, 
                      cuDoubleComplex *A,  magma_int_t lda, 
                      magma_int_t *info);
void magmablas_zlat2c(char uplo, magma_int_t n, 
                      cuDoubleComplex *A,  magma_int_t lda, 
                      cuFloatComplex  *SA, magma_int_t ldsa, 
                      magma_int_t *info);

#ifdef __cplusplus
}
#endif

#endif
