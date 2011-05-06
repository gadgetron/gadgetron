/*
 *   -- MAGMA (version 1.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2010
 *
 * @generated ds
 */

#ifndef _MAGMA_DS_H_
#define _MAGMA_DS_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_int_t magma_dsgesv_gpu(  char trans, magma_int_t N, magma_int_t NRHS, 
			       double *dA, magma_int_t ldda, 
			       magma_int_t *IPIV, magma_int_t *dIPIV,
			       double *dB, magma_int_t lddb, 
			       double *dX, magma_int_t lddx, 
			       double *dworkd, float *dworks,
			       magma_int_t *iter, magma_int_t *info);

magma_int_t magma_dsgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
			       float  *dA, magma_int_t ldda,
                               magma_int_t *ipiv, 
			       double *dB, magma_int_t lddb,
			       double *dX, magma_int_t lddx,
                               float  *dSX, 
			       magma_int_t *info );

magma_int_t magma_dsposv_gpu( char uplo, magma_int_t n, magma_int_t nrhs, 
                              double *dA, magma_int_t ldda, 
                              double *dB, magma_int_t lddb, 
                              double *dX, magma_int_t lddx, 
                              double *dworkd, float *dworks,
                              magma_int_t *iter, magma_int_t *info);

magma_int_t magma_dsgeqrsv_gpu(magma_int_t M, magma_int_t N, magma_int_t NRHS, 
			       double *dA,  magma_int_t ldda, 
			       double *dB,  magma_int_t lddb, 
			       double *dX,  magma_int_t lddx,
			       magma_int_t *iter,    magma_int_t *info);
  

#ifdef __cplusplus
}
#endif

#endif /* _MAGMA_D_H_ */
