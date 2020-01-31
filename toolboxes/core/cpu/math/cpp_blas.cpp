#include "cpp_blas.h"

#ifdef USE_MKL
#include "mkl.h"
#else
extern "C" {
#include "cblas.h"
}
#endif // MKL_FOUND

#ifdef OPENBLAS_SEQUENTIAL // Check for OpenBlas.
#define CBLAS_COMPLEX_FLOAT openblas_complex_float
#define CBLAS_COMPLEX_DOUBLE openblas_complex_double
#else
#define CBLAS_COMPLEX_FLOAT void
#define CBLAS_COMPLEX_DOUBLE void
#endif

float Gadgetron::BLAS::asum(size_t N, const float* x, size_t incx) {
    return cblas_sasum(N, x, incx);
}

double Gadgetron::BLAS::asum(size_t N, const double* x, size_t incx) {
    return cblas_dasum(N, x, incx);
}

float Gadgetron::BLAS::asum(size_t N, const std::complex<float>* x, size_t incx) {
    return cblas_scasum(N, (float*)x, incx);
}

double Gadgetron::BLAS::asum(size_t N, const std::complex<double>* x, size_t incx) {
    return cblas_dzasum(N, (double*)x, incx);
}

float Gadgetron::BLAS::asum(size_t N, const Gadgetron::complext<float>* x, size_t incx) {
    return cblas_scasum(N, (float*)x, incx);
}

double Gadgetron::BLAS::asum(size_t N, const Gadgetron::complext<double>* x, size_t incx) {
    return cblas_dzasum(N, (double*)x, incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const float* x, size_t incx) {
    return cblas_isamax(N, x, incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const double* x, size_t incx) {
    return cblas_idamax(N, x, incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const std::complex<float>* x, size_t incx) {
    return cblas_icamax(N, (float*)x, incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const std::complex<double>* x, size_t incx) {
    return cblas_izamax(N, (double*)x, incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const Gadgetron::complext<float>* x, size_t incx) {
    return cblas_icamax(N, (float*)x, incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const Gadgetron::complext<double>* x, size_t incx) {
    return cblas_izamax(N, (double*)x, incx);
}

void Gadgetron::BLAS::axpy(size_t N, float a, const float* x, size_t incx, float* y, size_t incy) {
    cblas_saxpy(N, a, x, incx, y, incy);
}

void Gadgetron::BLAS::axpy(size_t N, double a, const double* x, size_t incx, double* y, size_t incy) {
    cblas_daxpy(N, a, x, incx, y, incy);
}

void Gadgetron::BLAS::axpy(
    size_t N, std::complex<float> a, const std::complex<float>* x, size_t incx, std::complex<float>* y, size_t incy) {
    cblas_caxpy(N, (float*)&a, (float*)x, incx, (float*)y, incy);
}

void Gadgetron::BLAS::axpy(size_t N, std::complex<double> a, const std::complex<double>* x, size_t incx,
    std::complex<double>* y, size_t incy) {
    cblas_zaxpy(N, (double*)&a, (double*)x, incx, (double*)y, incy);
}

void Gadgetron::BLAS::axpy(size_t N, Gadgetron::complext<float> a, const Gadgetron::complext<float>* x, size_t incx,
    Gadgetron::complext<float>* y, size_t incy) {
    cblas_caxpy(N, (float*)&a, (float*)x, incx, (float*)y, incy);
}

void Gadgetron::BLAS::axpy(size_t N, Gadgetron::complext<double> a, const Gadgetron::complext<double>* x, size_t incx,
    Gadgetron::complext<double>* y, size_t incy) {
    cblas_zaxpy(N, (double*)&a, (double*)x, incx, (double*)y, incy);
}

float Gadgetron::BLAS::dot(size_t N, const float* x, size_t incx, const float* y, size_t incy) {
    return cblas_sdot(N, x, incx, y, incy);
}

double Gadgetron::BLAS::dot(size_t N, const double* x, size_t incx, const double* y, size_t incy) {
    return cblas_ddot(N, x, incx, y, incy);
}

std::complex<double> Gadgetron::BLAS::dot(
    size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy) {
    std::complex<double> result;
    cblas_zdotc_sub(N, (double*)x, incx, (double*)y, incy, (CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

std::complex<float> Gadgetron::BLAS::dot(
    size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy) {
    std::complex<float> result;
    cblas_cdotc_sub(N, (float*)x, incx, (float*)y, incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<float> Gadgetron::BLAS::dot(
    size_t N, const Gadgetron::complext<float>* x, size_t incx, const Gadgetron::complext<float>* y, size_t incy) {
    complext<float> result;
    cblas_cdotc_sub(N, (float*)x, incx, (float*)y, incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<double> Gadgetron::BLAS::dot(
    size_t N, const Gadgetron::complext<double>* x, size_t incx, const Gadgetron::complext<double>* y, size_t incy) {
    complext<double> result;
    cblas_zdotc_sub(N, (double*)x, incx, (double*)y, incy, (CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

std::complex<double> Gadgetron::BLAS::dotu(
    size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy) {
    std::complex<double> result;
    cblas_zdotu_sub(N, (double*)x, incx, (double*)y, incy, (CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

std::complex<float> Gadgetron::BLAS::dotu(
    size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy) {
    std::complex<float> result;
    cblas_cdotu_sub(N, (float*)x, incx, (float*)y, incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<float> Gadgetron::BLAS::dotu(
    size_t N, const Gadgetron::complext<float>* x, size_t incx, const Gadgetron::complext<float>* y, size_t incy) {
    complext<float> result;
    cblas_cdotu_sub(N, (float*)x, incx, (float*)y, incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<double> Gadgetron::BLAS::dotu(
    size_t N, const Gadgetron::complext<double>* x, size_t incx, const Gadgetron::complext<double>* y, size_t incy) {
    complext<double> result;
    cblas_zdotu_sub(N, (double*)x, incx, (double*)y, incy, (CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

std::complex<double> Gadgetron::BLAS::dotc(
    size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy) {
    std::complex<double> result;
    cblas_zdotc_sub(N, (double*)x, incx, (double*)y, incy, (CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

std::complex<float> Gadgetron::BLAS::dotc(
    size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy) {
    std::complex<float> result;
    cblas_cdotc_sub(N, (const float*)x, incx, (const float*)y, incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<float> Gadgetron::BLAS::dotc(
    size_t N, const Gadgetron::complext<float>* x, size_t incx, const Gadgetron::complext<float>* y, size_t incy) {
    complext<float> result;
    cblas_cdotc_sub(N, (const float*)x, incx, (const float*)y, incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<double> Gadgetron::BLAS::dotc(
    size_t N, const Gadgetron::complext<double>* x, size_t incx, const Gadgetron::complext<double>* y, size_t incy) {
    complext<double> result;
    cblas_zdotc_sub(N, (const double*)x, incx, (const double*)y, incy, (CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

float Gadgetron::BLAS::nrm2(size_t N, const float* x, size_t incx) {
    return cblas_snrm2(N, x, incx);
}

double Gadgetron::BLAS::nrm2(size_t N, const double* x, size_t incx) {
    return cblas_dnrm2(N, x, incx);
}

float Gadgetron::BLAS::nrm2(size_t N, const std::complex<float>* x, size_t incx) {
    return cblas_scnrm2(N, (float*)x, incx);
}

double Gadgetron::BLAS::nrm2(size_t N, const std::complex<double>* x, size_t incx) {
    return cblas_dznrm2(N, (double*)x, incx);
}

float Gadgetron::BLAS::nrm2(size_t N, const Gadgetron::complext<float>* x, size_t incx) {
    return cblas_scnrm2(N, (float*)x, incx);
}

double Gadgetron::BLAS::nrm2(size_t N, const Gadgetron::complext<double>* x, size_t incx) {
    return cblas_dznrm2(N, (double*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, float a, float* x, size_t incx) {
    cblas_sscal(N, a, x, incx);
}

void Gadgetron::BLAS::scal(size_t N, double a, double* x, size_t incx) {
    cblas_dscal(N, a, x, incx);
}

void Gadgetron::BLAS::scal(size_t N, std::complex<double> a, std::complex<double>* x, size_t incx) {
    cblas_zscal(N, (double*)&a, (double*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, std::complex<float> a, std::complex<float>* x, size_t incx) {
    cblas_cscal(N, (float*)&a, (float*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, Gadgetron::complext<float> a, Gadgetron::complext<float>* x, size_t incx) {
    cblas_cscal(N, (float*)&a, (float*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, Gadgetron::complext<double> a, Gadgetron::complext<double>* x, size_t incx) {
    cblas_zscal(N, (double*)&a, (double*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, double a, std::complex<double>* x, size_t incx) {
    cblas_zdscal(N, a, (double*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, float a, std::complex<float>* x, size_t incx) {
    cblas_csscal(N, a, (float*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, float a, Gadgetron::complext<float>* x, size_t incx) {
    cblas_csscal(N, a, (float*)x, incx);
}

void Gadgetron::BLAS::scal(size_t N, double a, Gadgetron::complext<double>* x, size_t incx) {
    cblas_zdscal(N, a, (double*)x, incx);
}
void Gadgetron::BLAS::gemm(bool transa, bool transb, size_t m, size_t n, size_t k, float alpha, const float* a,
    size_t lda, const float* b, size_t ldb, float beta, float* c, size_t ldc) {
    cblas_sgemm(CblasColMajor, transa ? CblasTrans : CblasNoTrans, transb ? CblasTrans : CblasNoTrans, m, n, k, alpha,
        a, lda, b, ldb, beta, c, ldc);
}
void Gadgetron::BLAS::gemm(bool transa, bool transb, size_t m, size_t n, size_t k, double alpha, const double* a,
    size_t lda, const double* b, size_t ldb, double beta, double* c, size_t ldc) {
    cblas_dgemm(CblasColMajor, transa ? CblasTrans : CblasNoTrans, transb ? CblasTrans : CblasNoTrans, m, n, k, alpha,
        a, lda, b, ldb, beta, c, ldc);
}
void Gadgetron::BLAS::gemm(bool transa, bool transb, size_t m, size_t n, size_t k, std::complex<float> alpha,
    const std::complex<float>* a, size_t lda, const std::complex<float>* b, size_t ldb, std::complex<float> beta,
    std::complex<float>* c, size_t ldc) {

    cblas_cgemm(CblasColMajor, transa ? CblasConjTrans : CblasNoTrans, transb ? CblasConjTrans : CblasNoTrans, m, n, k,
        &alpha, (float*)a, lda, (float*)b, ldb, &beta, (float*)c, ldc);
}
void Gadgetron::BLAS::gemm(bool transa, bool transb, size_t m, size_t n, size_t k, std::complex<double> alpha,
    const std::complex<double>* a, size_t lda, const std::complex<double>* b, size_t ldb, std::complex<double> beta,
    std::complex<double>* c, size_t ldc) {

    cblas_zgemm(CblasColMajor, transa ? CblasConjTrans : CblasNoTrans, transb ? CblasConjTrans : CblasNoTrans, m, n, k,
        &alpha, (double*)a, lda, (double*)b, ldb, &beta, (double*)c, ldc);
}

void Gadgetron::BLAS::gemm(bool transa, bool transb, size_t m, size_t n, size_t k, complext<float> alpha,
    const complext<float>* a, size_t lda, const complext<float>* b, size_t ldb, complext<float> beta,
    complext<float>* c, size_t ldc) {

    cblas_cgemm(CblasColMajor, transa ? CblasConjTrans : CblasNoTrans, transb ? CblasConjTrans : CblasNoTrans, m, n, k,
        &alpha, (float*)a, lda, (float*)b, ldb, &beta, (float*)c, ldc);
}
void Gadgetron::BLAS::gemm(bool transa, bool transb, size_t m, size_t n, size_t k, complext<double> alpha,
    const complext<double>* a, size_t lda, const complext<double>* b, size_t ldb, complext<double> beta,
    complext<double>* c, size_t ldc) {

    cblas_zgemm(CblasColMajor, transa ? CblasConjTrans : CblasNoTrans, transb ? CblasConjTrans : CblasNoTrans, m, n, k,
        &alpha, (double*)a, lda, (double*)b, ldb, &beta, (double*)c, ldc);
}
void Gadgetron::BLAS::syrk(bool upper, bool trans, size_t n, size_t k, float alpha, const float* a, size_t lda,
    float beta, float* c, size_t ldc) {
    cblas_ssyrk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasTrans : CblasNoTrans, n, k, alpha, a, lda,
        beta, c, ldc);
}
void Gadgetron::BLAS::syrk(bool upper, bool trans, size_t n, size_t k, double alpha, const double* a, size_t lda,
    double beta, double* c, size_t ldc) {
    cblas_dsyrk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasTrans : CblasNoTrans, n, k, alpha, a, lda,
        beta, c, ldc);
}
void Gadgetron::BLAS::syrk(bool upper, bool trans, size_t n, size_t k, std::complex<float> alpha,
    const std::complex<float>* a, size_t lda, std::complex<float> beta, std::complex<float>* c, size_t ldc) {

    cblas_csyrk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasConjTrans : CblasNoTrans, n, k, &alpha,
        (float*)a, lda, &beta, (float*)c, ldc);
}
void Gadgetron::BLAS::syrk(bool upper, bool trans, size_t n, size_t k, std::complex<double> alpha,
    const std::complex<double>* a, size_t lda, std::complex<double> beta, std::complex<double>* c, size_t ldc) {
    cblas_zsyrk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasConjTrans : CblasNoTrans, n, k, &alpha,
        (double*)a, lda, &beta, (double*)c, ldc);
}
void Gadgetron::BLAS::syrk(bool upper, bool trans, size_t n, size_t k, complext<float> alpha,
                           const complext<float>* a, size_t lda, complext<float> beta, complext<float>* c, size_t ldc) {

    cblas_csyrk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasConjTrans : CblasNoTrans, n, k, &alpha,
                (float*)a, lda, &beta, (float*)c, ldc);
}
void Gadgetron::BLAS::syrk(bool upper, bool trans, size_t n, size_t k, complext<double> alpha,
                           const complext<double>* a, size_t lda, complext<double> beta, complext<double>* c, size_t ldc) {
    cblas_zsyrk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasConjTrans : CblasNoTrans, n, k, &alpha,
                (double*)a, lda, &beta, (double*)c, ldc);
}
void Gadgetron::BLAS::herk(bool upper, bool trans, size_t n, size_t k, float alpha,
    const std::complex<float>* a, size_t lda, float beta, std::complex<float>* c, size_t ldc) {
    cblas_cherk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasConjTrans : CblasNoTrans, n, k, alpha,
                (float*)a, lda, beta, (float*)c, ldc);
}
void Gadgetron::BLAS::herk(bool upper, bool trans, size_t n, size_t k, double alpha,
    const std::complex<double>* a, size_t lda, double beta, std::complex<double>* c, size_t ldc) {
    cblas_zherk(CblasColMajor, upper ? CblasUpper : CblasLower, trans ? CblasConjTrans : CblasNoTrans, n, k, alpha,
                (double*)a, lda, beta, (double*)c, ldc);
}
