/**
 *
 * @file transpose.h
 *
 *  MAGMA (version 1.0) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  November 2010
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11
 *
 * Macro to transpose matrices before and after computation
 * in LU kernels
 *
 **/

#ifndef _MAGMA_TRANSPOSE_H_
#define _MAGMA_TRANSPOSE_H_

#define magmablas_sgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_sinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(float), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magmablas_stranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magmablas_sgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_sinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magmablas_stranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magmablas_dgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_dinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(double), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magmablas_dtranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magmablas_dgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_dinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magmablas_dtranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magmablas_cgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_cinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(cuFloatComplex), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magmablas_ctranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magmablas_cgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_cinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magmablas_ctranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magmablas_zgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_zinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(cuDoubleComplex), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magmablas_ztranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magmablas_zgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magmablas_zinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magmablas_ztranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#endif /* _MAGMA_TRANSPOSE_H_ */
