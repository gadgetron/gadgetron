/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated s

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_s
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasSgemm magmablas_sgemm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_slarfb_gpu( char side, char trans, char direct, char storev,
		  magma_int_t m, magma_int_t n, magma_int_t k,
		  float *dV,    magma_int_t ldv, 
		  float *dT,    magma_int_t ldt,
		  float *dC,    magma_int_t ldc, 
		  float *dwork, magma_int_t ldwork)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Univ. of California Berkeley
       November 2010

    Purpose
    =======
    ZLARFB applies a real block reflector H or its transpose H' to a
    REAL m by n matrix C, from the left.

    Arguments
    =========
    SIDE    (input) CHARACTER
            = 'L': apply H or H' from the Left
            = 'R': apply H or H' from the Right (Not implemented)

    TRANS   (input) CHARACTER
            = 'N': apply H  (No transpose)      (Not implemented)
            = 'C': apply H' (Conjugate transpose)

    DIRECT  (input) CHARACTER
            Indicates how H is formed from a product of elementary
            reflectors
            = 'F': H = H(1) H(2) . . . H(k) (Forward)
            = 'B': H = H(k) . . . H(2) H(1) (Backward)

    STOREV  (input) CHARACTER
            Indicates how the vectors which define the elementary
            reflectors are stored:
            = 'C': Columnwise
            = 'R': Rowwise

    M       (input) INTEGER
            The number of rows of the matrix C.

    N       (input) INTEGER
            The number of columns of the matrix C.

    K       (input) INTEGER
            The order of the matrix T (= the number of elementary
            reflectors whose product defines the block reflector).

    DV      (input) REAL array, dimension (LDV,K)
            The matrix V. See further details.

    LDV     (input) INTEGER
            The leading dimension of the array V. LDV >= max(1,M);

    DT      (input) REAL array, dimension (LDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    LDT     (input) INTEGER
            The leading dimension of the array T. LDT >= K.

    DC      (input/output) REAL array, dimension (LDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C.

    LDC     (input) INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    WORK    (workspace) REAL array, dimension (LDWORK,K)

    LDWORK  (input) INTEGER
            The leading dimension of the array WORK. LDWORK >= max(1,N);
    ===================================================================      */

    float c_zero    = MAGMA_S_ZERO;
    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;

    /* Function Body */
    if (m <= 0 || n <= 0) {
        return MAGMA_SUCCESS;
    }

    char transt;
    if (trans == 'N' || trans == 'n')
      transt = MagmaTrans;
    else
      transt = MagmaNoTrans;

    if ( ( side  == 'r' || side  == 'R') ) 
      {
        magma_int_t info = -1;
        fprintf(stderr, "The case (side == right) is not implemented\n");
        magma_xerbla(__func__, &info);
        return MAGMA_ERR_ILLEGAL_VALUE;
      }

    if ( storev == 'c' || storev == 'C') {
        /*
          if (n==1 && m%32==0){
          // This is used when we have to apply H on only one vector
          magmablas_sgemvt(m, k, 1., dv_ref(0,0), ldv, dc_ref(0, 0), dwork);
          printf("m= %d, n = %d, ldwork = %d\n", m, k, ldwork);
          }
          else
        */
        cublasSgemm( MagmaTrans, MagmaNoTrans,
                     n, k, m,
                     c_one,  dC,    ldc,
                             dV,    ldv,
                     c_zero, dwork, ldwork);

        if (direct == 'F' || direct =='f')
            cublasStrmm( MagmaRight, MagmaUpper, transt, MagmaNonUnit,
                         n, k, 
                         c_one, dT,    ldt, 
                                dwork, ldwork);
        else
            cublasStrmm( MagmaRight, MagmaLower, transt, MagmaNonUnit,
                         n, k, 
                         c_one, dT,    ldt, 
                                dwork, ldwork);

        cublasSgemm( MagmaNoTrans, MagmaTrans, 
                     m, n, k, 
                     c_neg_one, dV,    ldv,
                                dwork, ldwork, 
                     c_one,     dC,    ldc);
    }
    else {
        cublasSgemm( MagmaNoTrans, MagmaTrans, 
                     m, k, n, 
                     c_one,  dC,    ldc,
                             dV,    ldv, 
                     c_zero, dwork, ldwork);

        cublasStrmm( MagmaRight, MagmaUpper, transt, MagmaNonUnit,
                     m, k, 
                     c_one, dT,    ldt, 
                            dwork, ldwork);

        cublasSgemm( MagmaNoTrans, MagmaNoTrans, 
                     m, n, k, 
                     c_neg_one, dwork, ldwork,
                                dV,    ldv,
                     c_one,     dC,    ldc);
    }
    return MAGMA_SUCCESS;
} /* magma_slarfb */
