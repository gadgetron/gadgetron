//
// Created by dchansen on 1/31/20.
//

#pragma once
#include <complex>
#include <complext.h>

namespace Gadgetron {
    namespace Lapack {
        using Int = int64_t;

        Int potrf(bool upper, size_t n, float* a, size_t lda);
        Int potrf(bool upper, size_t n, double* a, size_t lda);
        Int potrf(bool upper, size_t n, std::complex<float>* a, size_t lda);
        Int potrf(bool upper, size_t n, std::complex<double>* a, size_t lda);


        Int heev(bool eigenvectors, bool upper, size_t n, std::complex<float>* a, size_t lda, float* w);
        Int heev(bool eigenvectors, bool upper, size_t n, std::complex<double>* a, size_t lda, double* w);
        Int syev(bool eigenvectors, bool upper, size_t n, float* a, size_t lda, float* w);
        Int syev(bool eigenvectors, bool upper, size_t n, double* a, size_t lda, double* w);

        Int potri(bool upper, size_t n, float* a, size_t lda);
        Int potri(bool upper, size_t n, double* a, size_t lda);
        Int potri(bool upper, size_t n, std::complex<float>* a, size_t lda);
        Int potri(bool upper, size_t n, std::complex<double>* a, size_t lda);

        Int tritri(bool upper, bool unittriangular, size_t n, float* a, size_t lda);
        Int tritri(bool upper, bool unittriangular, size_t n, double* a, size_t lda);
        Int tritri(bool upper, bool unittriangular, size_t n, std::complex<float>* a, size_t lda);
        Int tritri(bool upper, bool unittriangular, size_t n, std::complex<double>* a, size_t lda);

        Int posv(bool upper, size_t n, size_t nrhs, float* a, size_t lda, float* b, size_t ldb);
        Int posv(bool upper, size_t n, size_t nrhs, double* a, size_t lda, double* b, size_t ldb);
        Int posv(bool upper, size_t n, size_t nrhs, std::complex<float>* a, size_t lda, std::complex<float>* b, size_t ldb);
        Int posv(bool upper, size_t n, size_t nrhs, std::complex<double>* a, size_t lda, std::complex<double>* b, size_t ldb);

        Int hesv(bool upper, size_t n, size_t nrhs, std::complex<float>* a, size_t lda, Int* ipiv, std::complex<float>*b, size_t ldb);
        Int hesv(bool upper, size_t n, size_t nrhs, std::complex<double>* a, size_t lda, Int* ipiv, std::complex<double>*b, size_t ldb);


        Int sysv(bool upper, size_t n, size_t nrhs, float* a, size_t lda, Int* ipiv, float*b, size_t ldb);
        Int sysv(bool upper, size_t n, size_t nrhs, double* a, size_t lda, Int* ipiv, double*b, size_t ldb);
        Int syv(bool upper, size_t n, size_t nrhs, std::complex<float>* a, size_t lda, Int* ipiv, std::complex<float>*b, size_t ldb);
        Int sysv(bool upper, size_t n, size_t nrhs, std::complex<double>* a, size_t lda, Int* ipiv, std::complex<double>*b, size_t ldb);
    }
}