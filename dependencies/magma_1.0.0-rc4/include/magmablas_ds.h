/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated ds
 */

#ifndef _MAGMABLAS_DS_H_
#define _MAGMABLAS_DS_H_

#include "cublas.h"
#include "cuda.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magmablas_dsaxpycp(float *, double *, magma_int_t, double *, double *);
void magmablas_daxpycp(double *, double *, magma_int_t, double *);
void magmablas_dslaswp(magma_int_t, double *, magma_int_t, float *, magma_int_t, magma_int_t *, magma_int_t);
void magmablas_dlag2s(magma_int_t M, magma_int_t N, const double *A, magma_int_t lda,  float *SA, magma_int_t ldsa, magma_int_t *info);

void magmablas_slag2d(magma_int_t M, magma_int_t N, 
                      float  *SA, magma_int_t ldsa, 
                      double *A,  magma_int_t lda, 
                      magma_int_t *info);
void magmablas_dlat2s(char uplo, magma_int_t n, 
                      double *A,  magma_int_t lda, 
                      float  *SA, magma_int_t ldsa, 
                      magma_int_t *info);

#ifdef __cplusplus
}
#endif

#endif
