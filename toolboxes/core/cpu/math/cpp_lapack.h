//
// Created by dchansen on 1/31/20.
//

#pragma once
#include <complex>
#include <complext.h>

namespace Gadgetron {
    namespace Lapack {
        long long potrf(bool upper, size_t n, float* a, size_t lda);
        long long potrf(bool upper, size_t n, double* a, size_t lda);
        long long potrf(bool upper, size_t n, std::complex<float>* a, size_t lda);
        long long potrf(bool upper, size_t n, std::complex<double>* a, size_t lda);


        long long heev(bool eigenvectors, bool upper, size_t n, std::complex<float>* a, size_t lda, float* w);
        long long heev(bool eigenvectors, bool upper, size_t n, std::complex<double>* a, size_t lda, double* w);
        long long syev(bool eigenvectors, bool upper, size_t n, float* a, size_t lda, float* w);
        long long syev(bool eigenvectors, bool upper, size_t n, double* a, size_t lda, double* w);

        long long potri(bool upper, size_t n, float* a, size_t lda);
        long long potri(bool upper, size_t n, double* a, size_t lda);
        long long potri(bool upper, size_t n, std::complex<float>* a, size_t lda);
        long long potri(bool upper, size_t n, std::complex<double>* a, size_t lda);

        long long tritri(bool upper, bool unittriangular, size_t n, float* a, size_t lda);
        long long tritri(bool upper, bool unittriangular, size_t n, double* a, size_t lda);
        long long tritri(bool upper, bool unittriangular, size_t n, std::complex<float>* a, size_t lda);
        long long tritri(bool upper, bool unittriangular, size_t n, std::complex<double>* a, size_t lda);

        long long posv(bool upper, size_t n, size_t nrhs, float* a, size_t lda, float* b, size_t ldb);
        long long posv(bool upper, size_t n, size_t nrhs, double* a, size_t lda, double* b, size_t ldb);
        long long posv(bool upper, size_t n, size_t nrhs, std::complex<float>* a, size_t lda, std::complex<float>* b, size_t ldb);
        long long posv(bool upper, size_t n, size_t nrhs, std::complex<double>* a, size_t lda, std::complex<double>* b, size_t ldb);
    }
}