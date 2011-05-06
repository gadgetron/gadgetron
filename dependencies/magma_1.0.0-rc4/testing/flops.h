/**
 *
 * @file flops.h
 *
 *  File provided by Univ. of Tennessee,
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2010-12-20
 *
 **/
/*
 * This file provide the flops formula for all Level 3 BLAS and some
 * Lapack routines.  Each macro uses the same size parameters as the
 * function associated and provide one formula for additions and one
 * for multiplications. Ecample to use these macros:
 *  - In real:
 *    flops = FMULS_GEMM((double)m, (double(n), (double(k)) 
 *          + FADDS_GEMM((double)m, (double(n), (double(k));
 *
 *  - In complex:
 *    flops = 6.0 * FMULS_GEMM((double)m, (double(n), (double(k)) 
 *          + 2.0 * FADDS_GEMM((double)m, (double(n), (double(k));
 *
 * All the formula are reported in the LAPACK Lawn 41:
 *     http://www.netlib.org/lapack/lawns/lawn41.ps
 */
#ifndef _FLOPS_H_
#define _FLOPS_H_

/*
 * Level 2 BLAS 
 */  
#define FMULS_GEMV(m, n) ((m) * (n) + 2. * (m))
#define FADDS_GEMV(m, n) ((m) * (n)           )

#define FMULS_SYMV(n) ((n) * (n) + 2. * (n))
#define FADDS_SYMV(n) ((n) * (n)           )
#define FMULS_HEMV FMULS_SYMV
#define FADDS_HEMV FADDS_SYMV

/*
 * Level 3 BLAS 
 */
#define FMULS_GEMM(m, n, k) ((m) * (n) * (k))
#define FADDS_GEMM(m, n, k) ((m) * (n) * (k))

#define FMULS_SYMM_L(m, n) ((m) * (m) * (n))
#define FADDS_SYMM_L(m, n) ((m) * (m) * (n))
#define FMULS_HEMM_L FMULS_SYMM_L
#define FADDS_HEMM_L FADDS_SYMM_L

#define FMULS_SYMM_R(m, n) ((m) * (n) * (n))
#define FADDS_SYMM_R(m, n) ((m) * (n) * (n))
#define FMULS_HEMM_R FMULS_SYMM_R
#define FADDS_HEMM_R FADDS_SYMM_R

#define FMULS_SYRK(k, n) (0.5 * (k) * (n) * ((n)+1))
#define FADDS_SYRK(k, n) (0.5 * (k) * (n) * ((n)+1))
#define FMULS_HERK FMULS_SYRK
#define FADDS_HERK FADDS_SYRK

#define FMULS_SYR2K(k, n) ((k) * (n) * (n)      )
#define FADDS_SYR2K(k, n) ((k) * (n) * (n) + (n))
#define FMULS_HER2K FMULS_SYR2K
#define FADDS_HER2K FADDS_SYR2K

#define FMULS_TRMM_L(m, n) (0.5 * (n) * (m) * ((m)+1))
#define FADDS_TRMM_L(m, n) (0.5 * (n) * (m) * ((m)-1))

#define FMULS_TRMM_R(m, n) (0.5 * (m) * (n) * ((n)+1))
#define FADDS_TRMM_R(m, n) (0.5 * (m) * (n) * ((n)-1))

#define FMULS_TRSM_L(m, n) (0.5 * (n) * (m) * ((m)+1))
#define FADDS_TRSM_L(m, n) (0.5 * (n) * (m) * ((m)-1))

#define FMULS_TRSM_R(m, n) (0.5 * (m) * (n) * ((n)+1))
#define FADDS_TRSM_R(m, n) (0.5 * (m) * (n) * ((n)-1))

/*
 * Lapack
 */
#define FMULS_GETRF(m, n) ( ((m) < (n)) ? (0.5 * (m) * ((m) * ((n) - (1./3.) * (m) - 1. ) + (n)) + (2. / 3.) * (m)) \
 			    :             (0.5 * (n) * ((n) * ((m) - (1./3.) * (n) - 1. ) + (m)) + (2. / 3.) * (n)) )
#define FADDS_GETRF(m, n) ( ((m) < (n)) ? (0.5 * (m) * ((m) * ((n) - (1./3.) * (m)      ) - (n)) + (1. / 6.) * (m)) \
			    :             (0.5 * (n) * ((n) * ((m) - (1./3.) * (n)      ) - (m)) + (1. / 6.) * (n)) )

#define FMULS_GETRI(n) ( (n) * ((5. / 6.) + (n) * ((2. / 3.) * (n) + 0.5)) )
#define FADDS_GETRI(n) ( (n) * ((5. / 6.) + (n) * ((2. / 3.) * (n) - 1.5)) )

#define FMULS_GETRS(n, nrhs) ((nrhs) * (n) *  (n)      )
#define FADDS_GETRS(n, nrhs) ((nrhs) * (n) * ((n) - 1 ))

#define FMULS_POTRF(n) ((n) * (((1. / 6.) * (n) + 0.5) * (n) + (1. / 3.)))
#define FADDS_POTRF(n) ((n) * (((1. / 6.) * (n)      ) * (n) - (1. / 6.)))

#define FMULS_POTRI(n) ( (n) * ((2. / 3.) + (n) * ((1. / 3.) * (n) + 1. )) )
#define FADDS_POTRI(n) ( (n) * ((1. / 6.) + (n) * ((1. / 3.) * (n) - 0.5)) )

#define FMULS_POTRS(n, nrhs) ((nrhs) * (n) * ((n) + 1 ))
#define FADDS_POTRS(n, nrhs) ((nrhs) * (n) * ((n) - 1 ))

//SPBTRF
//SPBTRS
//SSYTRF
//SSYTRI
//SSYTRS

#define FMULS_GEQRF(m, n) (((m) > (n)) ? ((n) * ((n) * (  0.5-(1./3.) * (n) + (m)) +    (m) + 23. / 6.)) \
                                       : ((m) * ((m) * ( -0.5-(1./3.) * (m) + (n)) + 2.*(n) + 23. / 6.)) )
#define FADDS_GEQRF(m, n) (((m) > (n)) ? ((n) * ((n) * (  0.5-(1./3.) * (n) + (m))          +  5. / 6.)) \
                                       : ((m) * ((m) * ( -0.5-(1./3.) * (m) + (n)) +    (n) +  5. / 6.)) )

#define FMULS_GEQLF(m, n) FMULS_GEQRF(m, n)
#define FADDS_GEQLF(m, n) FADDS_GEQRF(m, n)

#define FMULS_GERQF(m, n) (((m) > (n)) ? ((n) * ((n) * (  0.5-(1./3.) * (n) + (m)) +    (m) + 29. / 6.)) \
                                       : ((m) * ((m) * ( -0.5-(1./3.) * (m) + (n)) + 2.*(n) + 29. / 6.)) )
#define FADDS_GERQF(m, n) (((m) > (n)) ? ((n) * ((n) * ( -0.5-(1./3.) * (n) + (m)) +    (m) +  5. / 6.)) \
                                       : ((m) * ((m) * (  0.5-(1./3.) * (m) + (n)) +        +  5. / 6.)) )

#define FMULS_GELQF(m, n) FMULS_GERQF(m, n)
#define FADDS_GELQF(m, n) FADDS_GERQF(m, n)

#define FMULS_UNGQR(m, n, k) ((k) * (2.* (m) * (n) +  2. * (n) - 5./3. + (k) * ( 2./3. * (k) - ((m) + (n)) - 1.)))
#define FADDS_UNGQR(m, n, k) ((k) * (2.* (m) * (n) + (n) - (m) + 1./3. + (k) * ( 2./3. * (k) - ((m) + (n))     )))
#define FMULS_UNGQL FMULS_UNGQR
#define FMULS_ORGQR FMULS_UNGQR
#define FMULS_ORGQL FMULS_UNGQR
#define FADDS_UNGQL FADDS_UNGQR
#define FADDS_ORGQR FADDS_UNGQR
#define FADDS_ORGQL FADDS_UNGQR

#define FMULS_UNGRQ(m, n, k) ((k) * (2.* (m) * (n) + (m) + (n) - 2./3. + (k) * ( 2./3. * (k) - ((m) + (n)) - 1.)))
#define FADDS_UNGRQ(m, n, k) ((k) * (2.* (m) * (n) + (m) - (n) + 1./3. + (k) * ( 2./3. * (k) - ((m) + (n))     )))
#define FMULS_UNGLQ FMULS_UNGRQ
#define FMULS_ORGRQ FMULS_UNGRQ
#define FMULS_ORGLQ FMULS_UNGRQ
#define FADDS_UNGLQ FADDS_UNGRQ
#define FADDS_ORGRQ FADDS_UNGRQ
#define FADDS_ORGLQ FADDS_UNGRQ

#define FMULS_GEQRS(m, n, nrhs) ((nrhs) * ((n) * ( 2.* (m) - 0.5 * (n) + 2.5)))
#define FADDS_GEQRS(m, n, nrhs) ((nrhs) * ((n) * ( 2.* (m) - 0.5 * (n) + 0.5)))

//UNMQR, UNMLQ, UNMQL, UNMRQ (Left)
//UNMQR, UNMLQ, UNMQL, UNMRQ (Right)

#define FMULS_TRTRI(n) ((n) * ((n) * ( 1./6. * (n) + 0.5 ) + 1./3.))
#define FADDS_TRTRI(n) ((n) * ((n) * ( 1./6. * (n) - 0.5 ) + 1./3.))

#define FMULS_GEHRD(n) ( (n) * ((n) * (5./3. *(n) + 0.5) - 7./6.) - 13. )
#define FADDS_GEHRD(n) ( (n) * ((n) * (5./3. *(n) - 1. ) - 2./3.) -  8. )

#define FMULS_SYTRD(n) ( (n) *  ( (n) * ( 2./3. * (n) + 2.5 ) - 1./6. ) - 15.)
#define FADDS_SYTRD(n) ( (n) *  ( (n) * ( 2./3. * (n) + 1.  ) - 8./3. ) -  4.)
#define FMULS_HETRD FMULS_SYTRD
#define FADDS_HETRD FADDS_SYTRD

#define FMULS_GEBRD(m, n) ( ((m) >= (n)) ? ((n) * ((n) * (2. * (m) - 2./3. * (n) + 2. )       + 20./3.)) \
                            :              ((m) * ((m) * (2. * (n) - 2./3. * (m) + 2. )       + 20./3.)) )
#define FADDS_GEBRD(m, n) ( ((m) >= (n)) ? ((n) * ((n) * (2. * (m) - 2./3. * (n) + 1. ) - (m) +  5./3.)) \
                            :              ((m) * ((m) * (2. * (n) - 2./3. * (m) + 1. ) - (n) +  5./3.)) )

#endif /* _FLOPS_H_ */
