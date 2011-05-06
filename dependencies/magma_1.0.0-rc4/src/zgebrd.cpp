/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/
#include "common_magma.h"

#define A(i, j)  (a + (j)*lda  + (i))
#define dA(i, j) (da+ (j)*ldda + (i))

// === Define what BLAS to use ============================================
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasZgemm magmablas_zgemm
#endif
// === End defining what BLAS to use ======================================

extern "C" magma_int_t
magma_zgebrd(magma_int_t m, magma_int_t n,
             cuDoubleComplex *a, magma_int_t lda, double *d, double *e,
             cuDoubleComplex *tauq, cuDoubleComplex *taup, 
             cuDoubleComplex *work, magma_int_t lwork, 
	     magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    ZGEBRD reduces a general complex M-by-N matrix A to upper or lower
    bidiagonal form B by an orthogonal transformation: Q\*\*H * A * P = B.

    If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows in the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns in the matrix A.  N >= 0.

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)
            On entry, the M-by-N general matrix to be reduced.
            On exit,
            if m >= n, the diagonal and the first superdiagonal are
              overwritten with the upper bidiagonal matrix B; the
              elements below the diagonal, with the array TAUQ, represent
              the orthogonal matrix Q as a product of elementary
              reflectors, and the elements above the first superdiagonal,
              with the array TAUP, represent the orthogonal matrix P as
              a product of elementary reflectors;
            if m < n, the diagonal and the first subdiagonal are
              overwritten with the lower bidiagonal matrix B; the
              elements below the first subdiagonal, with the array TAUQ,
              represent the orthogonal matrix Q as a product of
              elementary reflectors, and the elements above the diagonal,
              with the array TAUP, represent the orthogonal matrix P as
              a product of elementary reflectors.
            See Further Details.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    D       (output) COMPLEX_16 array, dimension (min(M,N))
            The diagonal elements of the bidiagonal matrix B:
            D(i) = A(i,i).

    E       (output) COMPLEX_16 array, dimension (min(M,N)-1)
            The off-diagonal elements of the bidiagonal matrix B:
            if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
            if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.

    TAUQ    (output) COMPLEX_16 array dimension (min(M,N))
            The scalar factors of the elementary reflectors which
            represent the orthogonal matrix Q. See Further Details.

    TAUP    (output) COMPLEX_16 array, dimension (min(M,N))
            The scalar factors of the elementary reflectors which
            represent the orthogonal matrix P. See Further Details.

    WORK    (workspace/output) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The length of the array WORK.  LWORK >= max(1,M,N).
            For optimum performance LWORK >= (M+N)*NB, where NB
            is the optimal blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value.

    Further Details
    ===============
    The matrices Q and P are represented as products of elementary
    reflectors:

    If m >= n,
       Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
    Each H(i) and G(i) has the form:
       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
    where tauq and taup are complex scalars, and v and u are complex vectors;
    v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
    u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
    tauq is stored in TAUQ(i) and taup in TAUP(i).

    If m < n,
       Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
    Each H(i) and G(i) has the form:
       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
    where tauq and taup are complex scalars, and v and u are complex vectors;
    v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
    u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
    tauq is stored in TAUQ(i) and taup in TAUP(i).

    The contents of A on exit are illustrated by the following examples:

    m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):

      (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
      (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
      (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
      (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
      (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
      (  v1  v2  v3  v4  v5 )

    where d and e denote diagonal and off-diagonal elements of B, vi
    denotes an element of the vector defining H(i), and ui an element of
    the vector defining G(i).
    =====================================================================    */


    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex *da, *dwork;

    magma_int_t ncol, nrow, jmax, nb, ldda;

    static magma_int_t i, j, nx;
    static cuDoubleComplex ws;
    static magma_int_t iinfo;

    magma_int_t minmn;
    magma_int_t ldwrkx, ldwrky, lwkopt;
    magma_int_t lquery;

    nb   = magma_get_zgebrd_nb(n);
    ldda = m;

    lwkopt = (m + n) * nb;
    work[0] = MAGMA_Z_MAKE( lwkopt, 0. );
    lquery = (lwork == -1);
    
    /* Check arguments */
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max(1,m)) {
	*info = -4;
    } else if ( (lwork < max( max(1, m), n)) && (! lquery) ) {
        *info = -10;
    }
    if (*info < 0)
	return MAGMA_ERR_ILLEGAL_VALUE;
    else if (lquery)
	return MAGMA_SUCCESS;

    /* Quick return if possible */
    minmn = min(m,n);
    if (minmn == 0) {
        work[0] = c_one;
        return MAGMA_SUCCESS;
    }

    if ( CUBLAS_STATUS_SUCCESS 
         != cublasAlloc(n*ldda+(m+n)*nb, sizeof(cuDoubleComplex), (void**)&da) ) {
        fprintf (stderr, "!!!! device memory allocation error in zgebrd\n" );
        return MAGMA_ERR_CUBLASALLOC; 
    }
    dwork = da + (n)*ldda;

    MAGMA_Z_SET2REAL( ws, max(m,n) );
    ldwrkx = m;
    ldwrky = n;

    /* Set the block/unblock crossover point NX. */
    nx = 128;

    /* Copy the matrix to the GPU */
    if (minmn-nx>=1)
        cublasSetMatrix(m, n, sizeof(cuDoubleComplex), a, lda, da, ldda);

    for (i=0; i< (minmn - nx); i += nb) {

        /*  Reduce rows and columns i:i+nb-1 to bidiagonal form and return
            the matrices X and Y which are needed to update the unreduced
            part of the matrix */
        nrow = m - i;
        ncol = n - i;

        /*   Get the current panel (no need for the 1st iteration) */
        if ( i > 0 ) {
            cublasGetMatrix(nrow, nb, sizeof(cuDoubleComplex),
                            dA(i, i), ldda, 
                            A( i, i), lda );
            cublasGetMatrix(nb, ncol - nb, sizeof(cuDoubleComplex),
                            dA(i, i+nb), ldda,
                            A( i, i+nb), lda);
        }

        magma_zlabrd_gpu(nrow, ncol, nb,
                         A(i, i),          lda,    dA(i, i),          ldda,
                         d+i, e+i, tauq+i, taup+i,
                         work,             ldwrkx, dwork,             ldwrkx,  // x, dx
                         work+(ldwrkx*nb), ldwrky, dwork+(ldwrkx*nb), ldwrky); // y, dy

        /*  Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
            of the form  A := A - V*Y' - X*U' */
        nrow = m - i - nb;
        ncol = n - i - nb;

        // Send Y back to the GPU
	cublasSetMatrix(nrow, nb, sizeof(cuDoubleComplex),
			work  + nb, ldwrkx,
			dwork + nb, ldwrkx);
	cublasSetMatrix(ncol, nb, sizeof(cuDoubleComplex),
			work  + (ldwrkx+1)*nb, ldwrky,
			dwork + (ldwrkx+1)*nb, ldwrky);

        cublasZgemm( MagmaNoTrans, MagmaConjTrans, 
                     nrow, ncol, nb, 
                     c_neg_one, dA(i+nb, i   ),      ldda,
                                dwork+(ldwrkx+1)*nb, ldwrky,
                     c_one,     dA(i+nb, i+nb),      ldda);

        cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                     nrow, ncol, nb, 
                     c_neg_one, dwork+nb,         ldwrkx,
                                dA( i,    i+nb ), ldda,
                     c_one,     dA( i+nb, i+nb ), ldda);

        /* Copy diagonal and off-diagonal elements of B back into A */
        if (m >= n) {
            jmax = i + nb;
            for (j = i; j < jmax; ++j) {
                *A(j, j  ) = MAGMA_Z_MAKE( d[j], 0. );
                *A(j, j+1) = MAGMA_Z_MAKE( e[j], 0. );
            }
        } else {
            jmax = i + nb;
            for (j = i; j < jmax; ++j) {
                *A(j,   j ) = MAGMA_Z_MAKE( d[j], 0. );
                *A(j+1, j ) = MAGMA_Z_MAKE( e[j], 0. );
                /* L20: */
            }
        }
    }

    /* Use unblocked code to reduce the remainder of the matrix */
    nrow = m - i;
    ncol = n - i;

    if ( 0 < (n-nx) )
        cublasGetMatrix(nrow, ncol, sizeof(cuDoubleComplex),
                        dA(i, i), ldda,
                        A( i, i), lda);

    lapackf77_zgebd2( &nrow, &ncol, 
                      A(i, i), &lda, d+i, e+i,
                      tauq+i, taup+i, work, &iinfo);
    work[0] = ws;

    cublasFree(da);
    return MAGMA_SUCCESS;
} /* zgebrd_ */

