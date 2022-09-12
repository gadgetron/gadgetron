#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"
#include <cpu/math/hoNDArray_linalg.h>
#include <gtest/gtest.h>
#include <random>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_blas_Real : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

typedef Types<float, double> realImplementations;

TYPED_TEST_SUITE(hoNDArray_blas_Real, realImplementations);

TYPED_TEST(hoNDArray_blas_Real,dotTest){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(),real(dot(&this->Array,&this->Array)));
  fill(&this->Array2,TypeParam(2));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements()*2,real(dot(&this->Array,&this->Array2)));
}

TYPED_TEST(hoNDArray_blas_Real,axpyTest){
  fill(&this->Array,TypeParam(71));
  fill(&this->Array2,TypeParam(97));
  axpy(TypeParam(11),&this->Array,&this->Array2);
  TypeParam val = this->Array2.get_data_ptr()[10];
  EXPECT_FLOAT_EQ(878,real(val));
}

TYPED_TEST(hoNDArray_blas_Real,nrm2Test){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(std::sqrt((double)this->Array.get_number_of_elements()),nrm2(&this->Array));
  fill(&this->Array,TypeParam(3));
  EXPECT_FLOAT_EQ(std::sqrt(3.0*3.0*this->Array.get_number_of_elements()),nrm2(&this->Array));
}

TYPED_TEST(hoNDArray_blas_Real,scal){
    fill(&this->Array,TypeParam(33));
    scal(2,this->Array);
    EXPECT_FLOAT_EQ(66,real(this->Array[2]));
}

TYPED_TEST(hoNDArray_blas_Real,asumTest){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(),real(asum(&this->Array)));
  fill(&this->Array,TypeParam(-3));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements()*3,real(asum(&this->Array)));
}

TYPED_TEST(hoNDArray_blas_Real,aminTest){
  fill(&this->Array,TypeParam(100));
  this->Array.get_data_ptr()[23]=TypeParam(-50);
  EXPECT_EQ(23,amin(&this->Array));
  this->Array.get_data_ptr()[48]=TypeParam(2);
  EXPECT_EQ(48,amin(&this->Array));
}

TYPED_TEST(hoNDArray_blas_Real,amaxTest){
  fill(&this->Array,TypeParam(1));
  this->Array.get_data_ptr()[23]=TypeParam(2);
  EXPECT_EQ(23,amax(&this->Array));
  this->Array.get_data_ptr()[48]=TypeParam(-50);
  EXPECT_EQ(48,amax(&this->Array));
}


template <typename T> class hoNDArray_blas_Cplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

typedef Types<std::complex<float>, std::complex<double>, float_complext, double_complext> cplxImplementations;

TYPED_TEST_SUITE(hoNDArray_blas_Cplx, cplxImplementations);

TYPED_TEST(hoNDArray_blas_Cplx,dotTest){
  fill(&this->Array,TypeParam(1,1));
  TypeParam res = dot(&this->Array,&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(1,-1)*TypeParam(1,1))*this->Array.get_number_of_elements(),real(res));
  EXPECT_FLOAT_EQ(0,imag(res));
  fill(&this->Array2,TypeParam(2,2));
  res = dot(&this->Array2,&this->Array2);
  EXPECT_FLOAT_EQ(real(TypeParam(2,-2)*TypeParam(2,2))*this->Array.get_number_of_elements(),real(res));
  EXPECT_FLOAT_EQ(0,imag(res));
  res = dot(&this->Array,&this->Array2);
  EXPECT_FLOAT_EQ(real(TypeParam(1,-1)*TypeParam(2,2))*this->Array.get_number_of_elements(),real(res));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,-1)*TypeParam(2,2))*this->Array.get_number_of_elements(),imag(res));
}

TYPED_TEST(hoNDArray_blas_Cplx,axpyTest){
  fill(&this->Array,TypeParam(71.1,23.3));
  fill(&this->Array2,TypeParam(97.9,654.2));
  axpy(TypeParam(11.4),&this->Array,&this->Array2);
  TypeParam got = this->Array2.get_data_ptr()[546];
  TypeParam wanted = TypeParam(71.1,23.3)*TypeParam(11.4)+TypeParam(97.9,654.2);
  EXPECT_FLOAT_EQ(real(wanted),real(got));
  EXPECT_FLOAT_EQ(imag(wanted),imag(got));
}

TYPED_TEST(hoNDArray_blas_Cplx,nrm2Test){
  fill(&this->Array,TypeParam(1,1));
  EXPECT_FLOAT_EQ(std::sqrt(real(TypeParam(1,-1)*TypeParam(1,1))*this->Array.get_number_of_elements()),nrm2(&this->Array));
  fill(&this->Array,TypeParam(3.24,7.4));
  // There will be rounding errors from the sum, so loosen comparison
  EXPECT_NEAR(std::sqrt(real(TypeParam(3.24,-7.4)*TypeParam(3.24,7.4))*this->Array.get_number_of_elements()),nrm2(&this->Array),0.001);
}

TYPED_TEST(hoNDArray_blas_Cplx,asumTest){
  fill(&this->Array,TypeParam(-3,1));
  EXPECT_NEAR(4*this->Array.get_number_of_elements(),asum(&this->Array),0.0001);
}

TYPED_TEST(hoNDArray_blas_Cplx,aminTest){
  fill(&this->Array,TypeParam(100,101));
  this->Array.get_data_ptr()[23]=TypeParam(-50,-51);
  EXPECT_EQ(23,amin(&this->Array));

  fill(&this->Array, TypeParam(100, 101));
  this->Array.get_data_ptr()[48]=TypeParam(2,100);
  EXPECT_EQ(48,amin(&this->Array));

  fill(&this->Array, TypeParam(100, 101));
  this->Array.get_data_ptr()[1000]=TypeParam(-2,-76);
  EXPECT_EQ(1000,amin(&this->Array));
}

TYPED_TEST(hoNDArray_blas_Cplx,amaxTest){
  fill(&this->Array,TypeParam(1,1));
  this->Array.get_data_ptr()[768]=TypeParam(4,4);
  EXPECT_EQ(768,amax(&this->Array));
  this->Array.get_data_ptr()[48]=TypeParam(6,1);
  EXPECT_EQ(768,amax(&this->Array));
  this->Array.get_data_ptr()[999]=TypeParam(-3,-6);
  EXPECT_EQ(999,amax(&this->Array));
}

TEST(BLASLevel3,gemm_squared){

    hoNDArray<std::complex<float>> A(127,113);
    std::uniform_real_distribution<float> dist;
    std::mt19937_64 engine{};

    for (auto& d : A) d = {dist(engine),dist(engine)};

    auto B = A;

    auto C = A;
    auto C2 = A;

    gemm(C,A,true,B,false);
    gemm(C2,A,true,A,false);

    ASSERT_EQ(C,C2);



}
