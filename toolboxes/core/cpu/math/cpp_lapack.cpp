//
// Created by dchansen on 1/31/20.
//

#include "cpp_lapack.h"



#ifdef USE_MKL
#include "mkl.h"
#else
extern "C" {
#include "lapacke.h"
}
#endif

long long Gadgetron::Lapack::potrf(bool upper, size_t n, float* a, size_t lda) {
    return LAPACKE_spotrf(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n,a,lda);
}
long long Gadgetron::Lapack::potrf(bool upper, size_t n, double* a, size_t lda) {
    return LAPACKE_dpotrf(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n,a,lda);
}
long long Gadgetron::Lapack::potrf(bool upper, size_t n, std::complex<float>* a, size_t lda) {
    return LAPACKE_cpotrf(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n, reinterpret_cast<lapack_complex_float*>(a),lda);
}
long long Gadgetron::Lapack::potrf(bool upper, size_t n, std::complex<double>* a, size_t lda) {
    return LAPACKE_zpotrf(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n, reinterpret_cast<lapack_complex_double*>(a),lda);
}

long long Gadgetron::Lapack::heev(
    bool eigenvectors, bool upper, size_t n, std::complex<float>* a, size_t lda, float* w) {
    return LAPACKE_cheev(LAPACK_COL_MAJOR,eigenvectors ? 'V' : 'N',upper ? 'U' : 'L',n, reinterpret_cast<lapack_complex_float*>(a),lda, w);
}
long long Gadgetron::Lapack::heev(
    bool eigenvectors, bool upper, size_t n, std::complex<double>* a, size_t lda, double* w) {

    return LAPACKE_zheev(LAPACK_COL_MAJOR,eigenvectors ? 'V' : 'N',upper ? 'U' : 'L',n, reinterpret_cast<lapack_complex_double*>(a),lda, w);
}
long long Gadgetron::Lapack::syev(bool eigenvectors, bool upper, size_t n, float* a, size_t lda, float* w) {
    return LAPACKE_ssyev(LAPACK_COL_MAJOR,eigenvectors ? 'V' : 'N',upper ? 'U' : 'L',n, a,lda, w);
}
long long Gadgetron::Lapack::syev(bool eigenvectors, bool upper, size_t n, double* a, size_t lda, double* w) {
    return LAPACKE_dsyev(LAPACK_COL_MAJOR,eigenvectors ? 'V' : 'N',upper ? 'U' : 'L',n, a,lda, w);
}
long long Gadgetron::Lapack::potri(bool upper, size_t n, float* a, size_t lda) {
    return LAPACKE_spotri(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n,a,lda);
}
long long Gadgetron::Lapack::potri(bool upper, size_t n, double* a, size_t lda) {
    return LAPACKE_dpotri(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n,a,lda);
}
long long Gadgetron::Lapack::potri(bool upper, size_t n, std::complex<float>* a, size_t lda) {
    return LAPACKE_cpotri(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n, reinterpret_cast<lapack_complex_float*>(a),lda);
}
long long Gadgetron::Lapack::potri(bool upper, size_t n, std::complex<double>* a, size_t lda) {
    return LAPACKE_zpotri(LAPACK_COL_MAJOR, upper ? 'U' : 'L',n, reinterpret_cast<lapack_complex_double*>(a),lda);
}
long long Gadgetron::Lapack::tritri(bool upper, bool unittriangular, size_t n, float* a, size_t lda) {
    return LAPACKE_strtri(LAPACK_COL_MAJOR,upper ? 'U' : 'L', unittriangular ? 'U' : 'N',n,a,lda);
}
long long Gadgetron::Lapack::tritri(bool upper, bool unittriangular, size_t n, double* a, size_t lda) {
    return LAPACKE_dtrtri(LAPACK_COL_MAJOR,upper ? 'U' : 'L', unittriangular ? 'U' : 'N',n,a,lda);
}
long long Gadgetron::Lapack::tritri(bool upper, bool unittriangular, size_t n, std::complex<float>* a, size_t lda) {
    return LAPACKE_ctrtri(LAPACK_COL_MAJOR,upper ? 'U' : 'L', unittriangular ? 'U' : 'N',n, reinterpret_cast<lapack_complex_float*>(a),lda);
}
long long Gadgetron::Lapack::tritri(bool upper, bool unittriangular, size_t n, std::complex<double>* a, size_t lda) {
    return LAPACKE_ztrtri(LAPACK_COL_MAJOR,upper ? 'U' : 'L', unittriangular ? 'U' : 'N',n, reinterpret_cast<lapack_complex_double*>(a),lda);
}
long long Gadgetron::Lapack::posv(bool upper, size_t n, size_t nrhs, float* a, size_t lda, float* b, size_t ldb) {
    return LAPACKE_sposv(LAPACK_COL_MAJOR,upper ? 'U' : 'L', n, nrhs,a,lda,b,ldb);
}
long long Gadgetron::Lapack::posv(bool upper, size_t n, size_t nrhs, double* a, size_t lda, double* b, size_t ldb) {
    return LAPACKE_dposv(LAPACK_COL_MAJOR,upper ? 'U' : 'L', n, nrhs,a,lda,b,ldb);
}
long long Gadgetron::Lapack::posv(
    bool upper, size_t n, size_t nrhs, std::complex<float>* a, size_t lda, std::complex<float>* b, size_t ldb) {
    return LAPACKE_cposv(LAPACK_COL_MAJOR,upper ? 'U' : 'L', n, nrhs, reinterpret_cast<lapack_complex_float*>(a),lda,
        reinterpret_cast<lapack_complex_float*>(b),ldb);
}
long long Gadgetron::Lapack::posv(
    bool upper, size_t n, size_t nrhs, std::complex<double>* a, size_t lda, std::complex<double>* b, size_t ldb) {
    return LAPACKE_zposv(LAPACK_COL_MAJOR,upper ? 'U' : 'L', n, nrhs, reinterpret_cast<lapack_complex_double*>(a),lda,
                         reinterpret_cast<lapack_complex_double*>(b),ldb);
}

long long Gadgetron::Lapack::hesv(bool upper, size_t n, size_t nrhs, std::complex<float> *a, size_t lda, size_t *ipiv,
                                  std::complex<float> *b, size_t ldb) {
    return LAPACKE_chesv(LAPACK_COL_MAJOR, upper ? 'U' : 'L', n, nrhs, reinterpret_cast<lapack_complex_float*>(a),lda, ipiv,
                         reinterpret_cast<lapack_complex_float*>(b),ldb);
}

long long Gadgetron::Lapack::hesv(bool upper, size_t n, size_t nrhs, std::complex<double> *a, size_t lda, size_t *ipiv,
                                  std::complex<double> *b, size_t ldb) {
    return LAPACKE_zhesv(LAPACK_COL_MAJOR, upper ? 'U' : 'L', n, nrhs, reinterpret_cast<lapack_complex_double*>(a),lda, ipiv,
                         reinterpret_cast<lapack_complex_double*>(b),ldb);
}

long long
Gadgetron::Lapack::sysv(bool upper, size_t n, size_t nrhs, float *a, size_t lda, size_t *ipiv, float *b, size_t ldb) {
    return LAPACKE_ssysv(LAPACK_COL_MAJOR,upper ? 'U' : 'L',n,nrhs,a,lda,ipiv,b,ldb);
}

long long
Gadgetron::Lapack::sysv(bool upper, size_t n, size_t nrhs, double *a, size_t lda, size_t *ipiv, double *b, size_t ldb) {
    return LAPACKE_dsysv(LAPACK_COL_MAJOR,upper ? 'U' : 'L',n,nrhs,a,lda,ipiv,b,ldb);
}

long long Gadgetron::Lapack::syv(bool upper, size_t n, size_t nrhs, std::complex<float> *a, size_t lda, size_t *ipiv,
                                 std::complex<float> *b, size_t ldb) {

    return LAPACKE_csysv(LAPACK_COL_MAJOR, upper ? 'U' : 'L', n, nrhs, reinterpret_cast<lapack_complex_float*>(a),lda, ipiv,
}

long long Gadgetron::Lapack::sysv(bool upper, size_t n, size_t nrhs, std::complex<double> *a, size_t lda, size_t *ipiv,
                                  std::complex<double> *b, size_t ldb) {

    return LAPACKE_zsysv(LAPACK_COL_MAJOR, upper ? 'U' : 'L', n, nrhs, reinterpret_cast<lapack_complex_float*>(a),lda, ipiv,
}
