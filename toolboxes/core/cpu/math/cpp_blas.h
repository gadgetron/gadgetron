#pragma once
#include "complext.h"
#include <complex>

namespace Gadgetron {
    namespace BLAS {

         float asum(size_t N, const float* x, size_t  incx);
         double asum(size_t N, const double* x, size_t  incx);
         float asum(size_t N, const std::complex<float>* x, size_t  incx);
         double asum(size_t N, const std::complex<double>* x, size_t  incx);
         float asum(size_t N, const complext<float>* x, size_t  incx);
         double asum(size_t N, const complext<double>* x, size_t  incx);


         size_t amax(size_t N, const float* x, size_t incx);
         size_t amax(size_t N, const double* x, size_t incx);
         size_t amax(size_t N, const std::complex<float>* x, size_t incx);
         size_t amax(size_t N, const std::complex<double>* x, size_t incx);
         size_t amax(size_t N, const complext<float>* x, size_t incx);
         size_t amax(size_t N, const complext<double>* x, size_t incx);


         void axpy(size_t N, float a, const float* x, size_t incx, float* y, size_t incy);
         void axpy(size_t N, double a, const double* x, size_t incx, double* y, size_t incy);
         void axpy(size_t N, std::complex<float> a, const std::complex<float>* x, size_t incx, std::complex<float>* y, size_t incy);
         void axpy(size_t N, std::complex<double> a, const std::complex<double>* x, size_t incx, std::complex<double>* y, size_t incy);
         void axpy(size_t N, complext<float> a, const complext<float>* x, size_t incx, complext<float>* y, size_t incy);
         void axpy(size_t N, complext<double> a, const complext<double>* x, size_t incx, complext<double>* y, size_t incy);

         float dot(size_t N, const float* x, size_t incx, const float* y, size_t incy);
         double dot(size_t N, const double* x, size_t incx, const double* y, size_t incy);
         std::complex<double> dot(size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy);
         std::complex<float> dot(size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy);
         complext<float> dot(size_t N, const complext<float>* x, size_t incx, const complext<float>* y, size_t incy);
         complext<double> dot(size_t N, const complext<double>* x, size_t incx, const complext<double>* y, size_t incy);

         std::complex<double> dotu(size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy);
         std::complex<float> dotu(size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy);
         complext<float> dotu(size_t N, const complext<float>* x, size_t incx, const complext<float>* y, size_t incy);
         complext<double> dotu(size_t N, const complext<double>* x, size_t incx, const complext<double>* y, size_t incy);

         std::complex<double> dotc(size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy);
         std::complex<float> dotc(size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy);
         complext<float> dotc(size_t N, const complext<float>* x, size_t incx, const complext<float>* y, size_t incy);
         complext<double> dotc(size_t N, const complext<double>* x, size_t incx, const complext<double>* y, size_t incy);

         float nrm2(size_t N, const float* x, size_t  incx);
         double nrm2(size_t N, const double* x, size_t  incx);
         float nrm2(size_t N, const std::complex<float>* x, size_t  incx);
         double nrm2(size_t N, const std::complex<double>* x, size_t  incx);
         float nrm2(size_t N, const complext<float>* x, size_t  incx);
         double nrm2(size_t N, const complext<double>* x, size_t  incx);

         void scal(size_t N, float a, float* x, size_t incx);
         void scal(size_t N, double a, double* x, size_t incx);
         void scal(size_t N, std::complex<double> a, std::complex<double>* x, size_t incx);
         void scal(size_t N, std::complex<float> a, std::complex<float>* x, size_t incx);
         void scal(size_t N, complext<float> a, complext<float>* x, size_t incx);
         void scal(size_t N, complext<double> a, complext<double>* x, size_t incx);

         void scal(size_t N, double a, std::complex<double>* x, size_t incx);
         void scal(size_t N, float a, std::complex<float>* x, size_t incx);
         void scal(size_t N, float a, complext<float>* x, size_t incx);
         void scal(size_t N, double a, complext<double>* x, size_t incx);
//Level 3 BLAS routines
         void gemm(bool transa, bool transb, size_t m, size_t n, size_t k, float alpha, const float* a,size_t lda, const float* b, size_t ldb, float beta, float*c, size_t ldc);
         void gemm(bool transa, bool transb, size_t m, size_t n, size_t k, double alpha, const double* a,size_t lda, const double* b, size_t ldb, double beta, double*c, size_t ldc);
         void gemm(bool transa, bool transb, size_t m, size_t n, size_t k, std::complex<float> alpha, const std::complex<float>* a,size_t lda, const std::complex<float>* b, size_t ldb, std::complex<float> beta, std::complex<float>*c, size_t ldc);
         void gemm(bool transa, bool transb, size_t m, size_t n, size_t k, std::complex<double> alpha, const std::complex<double>* a,size_t lda, const std::complex<double>* b, size_t ldb, std::complex<double> beta, std::complex<double>*c, size_t ldc);
         void gemm(bool transa, bool transb, size_t m, size_t n, size_t k, complext<float> alpha, const complext<float>* a,size_t lda, const complext<float>* b, size_t ldb, complext<float> beta, complext<float>*c, size_t ldc);
        void gemm(bool transa, bool transb, size_t m, size_t n, size_t k, complext<double> alpha, const complext<double>* a,size_t lda, const complext<double>* b, size_t ldb, complext<double> beta, complext<double>*c, size_t ldc);

         void syrk(bool upper, bool trans, size_t n, size_t k, float alpha, const float* a, size_t lda, float beta, float* c, size_t ldc);
         void syrk(bool upper, bool trans, size_t n, size_t k, double alpha, const double* a, size_t lda, double beta, double* c, size_t ldc);
         void syrk(bool upper, bool trans, size_t n, size_t k, std::complex<float> alpha, const std::complex<float>* a, size_t lda, std::complex<float> beta, std::complex<float>* c, size_t ldc);
         void syrk(bool upper, bool trans, size_t n, size_t k, std::complex<double> alpha, const std::complex<double>* a, size_t lda, std::complex<double> beta, std::complex<double>* c, size_t ldc);
         void syrk(bool upper, bool trans, size_t n, size_t k, complext<float> alpha, const complext<float>* a, size_t lda, complext<float> beta, complext<float>* c, size_t ldc);
         void syrk(bool upper, bool trans, size_t n, size_t k, complext<double> alpha, const complext<double>* a, size_t lda, complext<double> beta, complext<double>* c, size_t ldc);


         void herk(bool upper, bool trans, size_t n, size_t k, float alpha, const std::complex<float>* a, size_t lda, float beta, std::complex<float>* c, size_t ldc);
         void herk(bool upper, bool trans, size_t n, size_t k, double alpha, const std::complex<double>* a, size_t lda, double beta, std::complex<double>* c, size_t ldc);



    }
}