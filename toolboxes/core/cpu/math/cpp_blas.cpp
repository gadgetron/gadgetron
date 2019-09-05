#include "cpp_blas.h"

#ifdef USE_MKL
    #include "mkl.h"
#else
    extern "C" {
        #include "cblas.h"
    }
#endif // MKL_FOUND

#ifdef OPENBLAS_SEQUENTIAL //Check for OpenBlas.
#define CBLAS_COMPLEX_FLOAT openblas_complex_float
#define CBLAS_COMPLEX_DOUBLE openblas_complex_double
#else
#define CBLAS_COMPLEX_FLOAT void
#define CBLAS_COMPLEX_DOUBLE void
#endif

float Gadgetron::BLAS::asum(size_t N, const float *x, size_t incx) {
    return cblas_sasum(N,x,incx);
}

double Gadgetron::BLAS::asum(size_t N, const double *x, size_t incx) {
    return cblas_dasum(N,x,incx);
}

float Gadgetron::BLAS::asum(size_t N, const std::complex<float> *x, size_t incx) {
    return cblas_scasum(N,(float*)x,incx);
}

double Gadgetron::BLAS::asum(size_t N, const std::complex<double> *x, size_t incx) {
    return cblas_dzasum(N,(double*)x,incx);
}

float Gadgetron::BLAS::asum(size_t N, const Gadgetron::complext<float> *x, size_t incx) {
    return cblas_scasum(N,(float*)x,incx);
}

double Gadgetron::BLAS::asum(size_t N, const Gadgetron::complext<double> *x, size_t incx) {
    return cblas_dzasum(N,(double*)x,incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const float *x, size_t incx) {
    return cblas_isamax(N,x,incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const double *x, size_t incx) {
    return cblas_idamax(N,x,incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const std::complex<float> *x, size_t incx) {
    return cblas_icamax(N, (float*) x,incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const std::complex<double> *x, size_t incx) {
    return cblas_izamax(N,(double*)x,incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const Gadgetron::complext<float> *x, size_t incx) {
    return cblas_icamax(N,(float*) x,incx);
}

size_t Gadgetron::BLAS::amax(size_t N, const Gadgetron::complext<double> *x, size_t incx) {
    return cblas_izamax(N,(double*)x,incx);
}


void Gadgetron::BLAS::axpy(size_t N, float a, const float *x, size_t incx, float* y, size_t incy) {
    cblas_saxpy(N,a,x,incx,y,incy);

}

void Gadgetron::BLAS::axpy(size_t N, double a, const double *x, size_t incx, double* y, size_t incy) {
    cblas_daxpy(N,a,x,incx,y,incy);

}

void Gadgetron::BLAS::axpy(size_t N, std::complex<float> a, const std::complex<float> *x, size_t incx,
                           std::complex<float>* y, size_t incy) {
    cblas_caxpy(N,(float*)&a,(float*)x,incx,(float*)y,incy);

}

void Gadgetron::BLAS::axpy(size_t N, std::complex<double> a, const std::complex<double> *x, size_t incx,
                           std::complex<double>* y, size_t incy) {
    cblas_zaxpy(N,(double*)&a,(double*)x,incx,(double*)y,incy);

}

void Gadgetron::BLAS::axpy(size_t N, Gadgetron::complext<float> a, const Gadgetron::complext<float> *x, size_t incx,
                           Gadgetron::complext<float>* y, size_t incy) {
    cblas_caxpy(N,(float*)&a,(float*)x,incx,(float*)y,incy);

}

void Gadgetron::BLAS::axpy(size_t N, Gadgetron::complext<double> a, const Gadgetron::complext<double> *x, size_t incx,
                           Gadgetron::complext<double>* y, size_t incy) {
    cblas_zaxpy(N,(double*)&a,(double*)x,incx,(double*)y,incy);

}

float Gadgetron::BLAS::dot(size_t N, const float *x, size_t incx, const float *y, size_t incy) {
    return cblas_sdot(N,x,incx,y,incy);
}

double Gadgetron::BLAS::dot(size_t N, const double *x, size_t incx, const double *y, size_t incy) {
    return cblas_ddot(N,x,incx,y,incy);
}

std::complex<double>
Gadgetron::BLAS::dot(size_t N, const std::complex<double> *x, size_t incx, const std::complex<double> *y, size_t incy) {
    std::complex<double> result;
    cblas_zdotc_sub(N, (double*)x,incx,(double*)y,incy,(CBLAS_COMPLEX_DOUBLE *)&result);
    return result;

}

std::complex<float>
Gadgetron::BLAS::dot(size_t N, const std::complex<float> *x, size_t incx, const std::complex<float> *y, size_t incy) {
    std::complex<float> result;
    cblas_cdotc_sub(N,(float*) x,incx,(float*) y,incy,(CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<float>
Gadgetron::BLAS::dot(size_t N, const Gadgetron::complext<float> *x, size_t incx, const Gadgetron::complext<float> *y,
                     size_t incy) {
    complext<float> result;
    cblas_cdotc_sub(N,(float*) x,incx,(float*) y,incy,(CBLAS_COMPLEX_FLOAT*) &result);
    return result;
}

Gadgetron::complext<double>
Gadgetron::BLAS::dot(size_t N, const Gadgetron::complext<double> *x, size_t incx, const Gadgetron::complext<double> *y,
                     size_t incy) {
    complext<double> result;
    cblas_zdotc_sub(N,(double*)x,incx,(double*)y,incy,(CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

std::complex<double>
Gadgetron::BLAS::dotu(size_t N, const std::complex<double> *x, size_t incx, const std::complex<double> *y,
                      size_t incy) {
    std::complex<double> result;
    cblas_zdotu_sub(N,(double*)x,incx,(double*)y,incy,(CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}

std::complex<float>
Gadgetron::BLAS::dotu(size_t N, const std::complex<float> *x, size_t incx, const std::complex<float> *y, size_t incy) {
    std::complex<float> result;
    cblas_cdotu_sub(N,(float*) x,incx,(float*) y,incy,(CBLAS_COMPLEX_FLOAT*) &result);
    return result;
}

Gadgetron::complext<float>
Gadgetron::BLAS::dotu(size_t N, const Gadgetron::complext<float> *x, size_t incx, const Gadgetron::complext<float> *y,
                      size_t incy) {
    complext<float> result;
    cblas_cdotu_sub(N,(float*) x,incx,(float*) y,incy,(CBLAS_COMPLEX_FLOAT*) &result);
    return result;
}

Gadgetron::complext<double>
Gadgetron::BLAS::dotu(size_t N, const Gadgetron::complext<double> *x, size_t incx, const Gadgetron::complext<double> *y,
                      size_t incy) {
    complext<double> result;
    cblas_zdotu_sub(N,(double*)x,incx,(double*)y,incy,(CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}


std::complex<double>
Gadgetron::BLAS::dotc(size_t N, const std::complex<double> *x, size_t incx, const std::complex<double> *y, size_t incy) {
    std::complex<double> result;
    cblas_zdotc_sub(N,(double*)x,incx,(double*)y,incy,(CBLAS_COMPLEX_DOUBLE*)&result);
    return result;

}

std::complex<float>
Gadgetron::BLAS::dotc(size_t N, const std::complex<float> *x, size_t incx, const std::complex<float> *y, size_t incy) {
    std::complex<float> result;
    cblas_cdotc_sub(N,(const float*) x,incx,(const float*) y,incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<float>
Gadgetron::BLAS::dotc(size_t N, const Gadgetron::complext<float> *x, size_t incx, const Gadgetron::complext<float> *y,
                     size_t incy) {
    complext<float> result;
    cblas_cdotc_sub(N,(const float*) x,incx,(const float*) y,incy, (CBLAS_COMPLEX_FLOAT*)&result);
    return result;
}

Gadgetron::complext<double>
Gadgetron::BLAS::dotc(size_t N, const Gadgetron::complext<double> *x, size_t incx, const Gadgetron::complext<double> *y,
                     size_t incy) {
    complext<double> result;
    cblas_zdotc_sub(N,(const double*)x,incx,(const double*)y,incy, (CBLAS_COMPLEX_DOUBLE*)&result);
    return result;
}


float Gadgetron::BLAS::nrm2(size_t N, const float *x, size_t incx) {
    return cblas_snrm2(N,x,incx);
}

double Gadgetron::BLAS::nrm2(size_t N, const double *x, size_t incx) {
    return cblas_dnrm2(N,x,incx);
}

float Gadgetron::BLAS::nrm2(size_t N, const std::complex<float> *x, size_t incx) {
    return cblas_scnrm2(N,(float*) x,incx);
}

double Gadgetron::BLAS::nrm2(size_t N, const std::complex<double> *x, size_t incx) {
    return cblas_dznrm2(N,(double*)x,incx);
}

float Gadgetron::BLAS::nrm2(size_t N, const Gadgetron::complext<float> *x, size_t incx) {
    return cblas_scnrm2(N,(float*) x,incx);
}

double Gadgetron::BLAS::nrm2(size_t N, const Gadgetron::complext<double> *x, size_t incx) {
    return cblas_dznrm2(N,(double*)x,incx);
}

void Gadgetron::BLAS::scal(size_t N, float a, float *x, size_t incx) {
    cblas_sscal(N,a,x,incx);
}

void Gadgetron::BLAS::scal(size_t N, double a, double *x, size_t incx) {
    cblas_dscal(N,a,x,incx);
}

void Gadgetron::BLAS::scal(size_t N, std::complex<double> a, std::complex<double> *x, size_t incx) {
    cblas_zscal(N,(double*)&a,(double*)x,incx);

}

void Gadgetron::BLAS::scal(size_t N, std::complex<float> a, std::complex<float> *x, size_t incx) {
    cblas_cscal(N,(float*) &a,(float*) x,incx);
}

void Gadgetron::BLAS::scal(size_t N, Gadgetron::complext<float> a, Gadgetron::complext<float> *x, size_t incx) {
    cblas_cscal(N,(float*) &a,(float*) x,incx);
}

void Gadgetron::BLAS::scal(size_t N, Gadgetron::complext<double> a, Gadgetron::complext<double> *x, size_t incx) {
    cblas_zscal(N,(double*)&a,(double*)x,incx);

}

void Gadgetron::BLAS::scal(size_t N, double a, std::complex<double> *x, size_t incx) {
    cblas_zdscal(N,a,(double*)x,incx);

}

void Gadgetron::BLAS::scal(size_t N, float a, std::complex<float> *x, size_t incx) {
    cblas_csscal(N,a,(float*) x,incx);

}

void Gadgetron::BLAS::scal(size_t N, float a, Gadgetron::complext<float> *x, size_t incx) {
    cblas_csscal(N,a,(float*) x,incx);

}

void Gadgetron::BLAS::scal(size_t N, double a, Gadgetron::complext<double> *x, size_t incx) {
    cblas_zdscal(N,a,(double*)x,incx);

}
