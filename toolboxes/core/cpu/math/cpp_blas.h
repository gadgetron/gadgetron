#pragma once
#include "cpucore_math_export.h"
#include "complext.h"
#include <complex>

namespace Gadgetron {
    namespace BLAS {

        EXPORTCPUCOREMATH float asum(size_t N, const float* x, size_t  incx);
        EXPORTCPUCOREMATH double asum(size_t N, const double* x, size_t  incx);
        EXPORTCPUCOREMATH float asum(size_t N, const std::complex<float>* x, size_t  incx);
        EXPORTCPUCOREMATH double asum(size_t N, const std::complex<double>* x, size_t  incx);
        EXPORTCPUCOREMATH float asum(size_t N, const complext<float>* x, size_t  incx);
        EXPORTCPUCOREMATH double asum(size_t N, const complext<double>* x, size_t  incx);


        EXPORTCPUCOREMATH size_t amax(size_t N, const float* x, size_t incx);
        EXPORTCPUCOREMATH size_t amax(size_t N, const double* x, size_t incx);
        EXPORTCPUCOREMATH size_t amax(size_t N, const std::complex<float>* x, size_t incx);
        EXPORTCPUCOREMATH size_t amax(size_t N, const std::complex<double>* x, size_t incx);
        EXPORTCPUCOREMATH size_t amax(size_t N, const complext<float>* x, size_t incx);
        EXPORTCPUCOREMATH size_t amax(size_t N, const complext<double>* x, size_t incx);


        EXPORTCPUCOREMATH void axpy(size_t N, float a, const float* x, size_t incx, float* y, size_t incy);
        EXPORTCPUCOREMATH void axpy(size_t N, double a, const double* x, size_t incx, double* y, size_t incy);
        EXPORTCPUCOREMATH void axpy(size_t N, std::complex<float> a, const std::complex<float>* x, size_t incx, std::complex<float>* y, size_t incy);
        EXPORTCPUCOREMATH void axpy(size_t N, std::complex<double> a, const std::complex<double>* x, size_t incx, std::complex<double>* y, size_t incy);
        EXPORTCPUCOREMATH void axpy(size_t N, complext<float> a, const complext<float>* x, size_t incx, complext<float>* y, size_t incy);
        EXPORTCPUCOREMATH void axpy(size_t N, complext<double> a, const complext<double>* x, size_t incx, complext<double>* y, size_t incy);

        EXPORTCPUCOREMATH float dot(size_t N, const float* x, size_t incx, const float* y, size_t incy);
        EXPORTCPUCOREMATH double dot(size_t N, const double* x, size_t incx, const double* y, size_t incy);
        EXPORTCPUCOREMATH std::complex<double> dot(size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy);
        EXPORTCPUCOREMATH std::complex<float> dot(size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy);
        EXPORTCPUCOREMATH complext<float> dot(size_t N, const complext<float>* x, size_t incx, const complext<float>* y, size_t incy);
        EXPORTCPUCOREMATH complext<double> dot(size_t N, const complext<double>* x, size_t incx, const complext<double>* y, size_t incy);

        EXPORTCPUCOREMATH std::complex<double> dotu(size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy);
        EXPORTCPUCOREMATH std::complex<float> dotu(size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy);
        EXPORTCPUCOREMATH complext<float> dotu(size_t N, const complext<float>* x, size_t incx, const complext<float>* y, size_t incy);
        EXPORTCPUCOREMATH complext<double> dotu(size_t N, const complext<double>* x, size_t incx, const complext<double>* y, size_t incy);

        EXPORTCPUCOREMATH std::complex<double> dotc(size_t N, const std::complex<double>* x, size_t incx, const std::complex<double>* y, size_t incy);
        EXPORTCPUCOREMATH std::complex<float> dotc(size_t N, const std::complex<float>* x, size_t incx, const std::complex<float>* y, size_t incy);
        EXPORTCPUCOREMATH complext<float> dotc(size_t N, const complext<float>* x, size_t incx, const complext<float>* y, size_t incy);
        EXPORTCPUCOREMATH complext<double> dotc(size_t N, const complext<double>* x, size_t incx, const complext<double>* y, size_t incy);

        EXPORTCPUCOREMATH float nrm2(size_t N, const float* x, size_t  incx);
        EXPORTCPUCOREMATH double nrm2(size_t N, const double* x, size_t  incx);
        EXPORTCPUCOREMATH float nrm2(size_t N, const std::complex<float>* x, size_t  incx);
        EXPORTCPUCOREMATH double nrm2(size_t N, const std::complex<double>* x, size_t  incx);
        EXPORTCPUCOREMATH float nrm2(size_t N, const complext<float>* x, size_t  incx);
        EXPORTCPUCOREMATH double nrm2(size_t N, const complext<double>* x, size_t  incx);

        EXPORTCPUCOREMATH void scal(size_t N, float a, float* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, double a, double* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, std::complex<double> a, std::complex<double>* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, std::complex<float> a, std::complex<float>* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, complext<float> a, complext<float>* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, complext<double> a, complext<double>* x, size_t incx);

        EXPORTCPUCOREMATH void scal(size_t N, double a, std::complex<double>* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, float a, std::complex<float>* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, float a, complext<float>* x, size_t incx);
        EXPORTCPUCOREMATH void scal(size_t N, double a, complext<double>* x, size_t incx);
    }
}