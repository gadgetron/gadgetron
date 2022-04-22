#include "cuNDArray_blas.h"
#include "cuNDArray_elemwise.h"

#include <gtest/gtest.h>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class cuNDArray_blas_Real : public ::testing::Test {
  protected:
    virtual void SetUp() {

        dims = std::vector<size_t>{37, 49, 23, 19}; // Using prime numbers for setup because they are messy

        Array = cuNDArray<T>(&dims);
        Array2 = cuNDArray<T>(&dims);
    }
    std::vector<size_t> dims;
    
    cuNDArray<T> Array;
    cuNDArray<T> Array2;
};

typedef Types<float, double> realImplementations;

TYPED_TEST_CASE(cuNDArray_blas_Real, realImplementations);

TYPED_TEST(cuNDArray_blas_Real, dotTest) {
    fill(&this->Array, TypeParam(1));
    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(), real(dot(&this->Array, &this->Array)));
    fill(&this->Array2, TypeParam(2));
    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements() * 2, real(dot(&this->Array, &this->Array2)));

    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements() * 2, real(dot(&this->Array, &this->Array2, 16)));
    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements() * 2, real(dot(&this->Array, &this->Array2, 17)));


}

TYPED_TEST(cuNDArray_blas_Real, axpyTest) {
    fill(&this->Array, TypeParam(71));
    fill(&this->Array2, TypeParam(97));
    axpy(TypeParam(11), &this->Array, &this->Array2);
    
    TypeParam val = this->Array2[10];
    EXPECT_FLOAT_EQ(878, real(val));

    fill(&this->Array2, TypeParam(97));
    axpy(TypeParam(11), &this->Array, &this->Array2 ,16);
    val = this->Array2[10];
    EXPECT_FLOAT_EQ(878, real(val));

    fill(&this->Array2, TypeParam(97));
    axpy(TypeParam(11), &this->Array, &this->Array2 ,17);
    val = this->Array2[10];
    EXPECT_FLOAT_EQ(878, real(val));

}

TYPED_TEST(cuNDArray_blas_Real, nrm2Test) {
    fill(&this->Array, TypeParam(1));
    EXPECT_FLOAT_EQ(std::sqrt((double)this->Array.get_number_of_elements()), nrm2(&this->Array));
    fill(&this->Array, TypeParam(3));
    EXPECT_FLOAT_EQ(std::sqrt(3.0 * 3.0 * this->Array.get_number_of_elements()), nrm2(&this->Array));

    // Some errors from the sum
    EXPECT_NEAR(std::sqrt(3.0 * 3.0 * this->Array.get_number_of_elements()), nrm2(&this->Array, 16),0.1);
    EXPECT_NEAR(std::sqrt(3.0 * 3.0 * this->Array.get_number_of_elements()), nrm2(&this->Array, 17),0.1);

}

TYPED_TEST(cuNDArray_blas_Real, asumTest) {
    fill(&this->Array, TypeParam(1));
    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(), real(asum(&this->Array)));
    fill(&this->Array, TypeParam(-3));
    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements() * 3, real(asum(&this->Array)));
    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements() * 3, real(asum(&this->Array,16)));
    EXPECT_FLOAT_EQ(this->Array.get_number_of_elements() * 3, real(asum(&this->Array,17)));

}

TYPED_TEST(cuNDArray_blas_Real, aminTest) {
    fill(&this->Array, TypeParam(100));
    TypeParam tmp(-50);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[4], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(4, amin(&this->Array));
    EXPECT_EQ(4, amin(&this->Array,16));
    EXPECT_EQ(4, amin(&this->Array,17));

    tmp = TypeParam(2);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[48], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(48, amin(&this->Array));
    EXPECT_EQ(48, amin(&this->Array,16));
    EXPECT_EQ(48, amin(&this->Array,17));



}

TYPED_TEST(cuNDArray_blas_Real, amaxTest) {
    fill(&this->Array, TypeParam(1));
    TypeParam tmp(2);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[4], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(4, amax(&this->Array));
    EXPECT_EQ(4, amax(&this->Array,16));
    EXPECT_EQ(4, amax(&this->Array,17));

    tmp = TypeParam(50);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[48], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(48, amax(&this->Array));
    EXPECT_EQ(48, amax(&this->Array,16));
    EXPECT_EQ(48, amax(&this->Array,17));


}

template <typename T> class cuNDArray_blas_Cplx : public ::testing::Test {
  protected:
    virtual void SetUp() {

        dims = std::vector<size_t>{37, 49}; // Using prime numbers for setup because they are messy
        Array = cuNDArray<T>(&dims);
        Array2 = cuNDArray<T>(&dims);

    }
    std::vector<size_t> dims;

    cuNDArray<T> Array;
    cuNDArray<T> Array2;

};

typedef Types</*std::complex<float>, std::complex<double>,*/ float_complext, double_complext> cplxImplementations;

TYPED_TEST_CASE(cuNDArray_blas_Cplx, cplxImplementations);

TYPED_TEST(cuNDArray_blas_Cplx, dotTest) {
    fill(&this->Array, TypeParam(1, 1));
    TypeParam res = dot(&this->Array, &this->Array);
    EXPECT_FLOAT_EQ(real(TypeParam(1, -1) * TypeParam(1, 1)) * this->Array.get_number_of_elements(), real(res));
    EXPECT_FLOAT_EQ(0, imag(res));
    fill(&this->Array2, TypeParam(2, 2));
    res = dot(&this->Array2, &this->Array2);
    EXPECT_FLOAT_EQ(real(TypeParam(2, -2) * TypeParam(2, 2)) * this->Array.get_number_of_elements(), real(res));
    EXPECT_FLOAT_EQ(0, imag(res));
    res = dot(&this->Array, &this->Array2);
    EXPECT_FLOAT_EQ(real(TypeParam(1, -1) * TypeParam(2, 2)) * this->Array.get_number_of_elements(), real(res));
    EXPECT_FLOAT_EQ(imag(TypeParam(1, -1) * TypeParam(2, 2)) * this->Array.get_number_of_elements(), imag(res));
}

TYPED_TEST(cuNDArray_blas_Cplx, axpyTest) {
    fill(&this->Array, TypeParam(71.1, 23.3));
    fill(&this->Array2, TypeParam(97.9, 654.2));
    axpy(TypeParam(11.4), &this->Array, &this->Array2);
    TypeParam got = this->Array2[546];
    TypeParam wanted = TypeParam(71.1, 23.3) * TypeParam(11.4) + TypeParam(97.9, 654.2);
    EXPECT_FLOAT_EQ(real(wanted), real(got));
    EXPECT_FLOAT_EQ(imag(wanted), imag(got));
}

TYPED_TEST(cuNDArray_blas_Cplx, nrm2Test) {
    fill(&this->Array, TypeParam(1, 1));
    EXPECT_FLOAT_EQ(std::sqrt(real(TypeParam(1, -1) * TypeParam(1, 1)) * this->Array.get_number_of_elements()),
                    nrm2(&this->Array));
    fill(&this->Array, TypeParam(3.24, 7.4));
    // There will be rounding errors from the sum, so loosen comparison
    EXPECT_NEAR(std::sqrt(real(TypeParam(3.24, -7.4) * TypeParam(3.24, 7.4)) * this->Array.get_number_of_elements()),
                nrm2(&this->Array), 0.001);
}

TYPED_TEST(cuNDArray_blas_Cplx, asumTest) {
    fill(&this->Array, TypeParam(-3, 1));
    EXPECT_NEAR(4 * this->Array.get_number_of_elements(), asum(&this->Array), 0.1);
}

TYPED_TEST(cuNDArray_blas_Cplx, aminTest) {
    fill(&this->Array, TypeParam(100, 101));
    TypeParam tmp(-50, -51);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[23], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(23, amin(&this->Array));
    EXPECT_EQ(23, amin(&this->Array,27));
    EXPECT_EQ(23, amin(&this->Array,17));

    tmp = TypeParam(2, 100);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[48], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(23, amin(&this->Array));
    EXPECT_EQ(23, amin(&this->Array,16));
    EXPECT_EQ(23, amin(&this->Array,17));
    tmp = TypeParam(-2, -76);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[1000], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(1000, amin(&this->Array));
    EXPECT_EQ(1000, amin(&this->Array,16));
    EXPECT_EQ(1000, amin(&this->Array,17));


}

TYPED_TEST(cuNDArray_blas_Cplx, amaxTest) {
    fill(&this->Array, TypeParam(1, 1));
    TypeParam tmp(4, 4);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[3], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(3, amax(&this->Array));
    EXPECT_EQ(3, amax(&this->Array,16));
    EXPECT_EQ(3, amax(&this->Array,17));

    tmp = TypeParam(6, 1);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[48], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(768, amax(&this->Array));
    EXPECT_EQ(768, amax(&this->Array,16));
    EXPECT_EQ(768, amax(&this->Array,17));

    tmp = TypeParam(-3, -6);
    CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[999], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));
    EXPECT_EQ(999, amax(&this->Array));
    EXPECT_EQ(999, amax(&this->Array,16));
    EXPECT_EQ(999, amax(&this->Array,17));

}
