/*
 * hoNDArray_blas_test.cpp
 *
 *  Created on: Feb 28, 2013
 *      Author: Dae
 */

#include "hoNDArray_blas.h"
#include "hoNDArray_elemwise.h"
#include <gtest/gtest.h>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_blas_Test : public ::testing::Test {
protected:
  virtual void SetUp() {
    unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
    Array = hoNDArray<T>(&dims);
    Array2 = hoNDArray<T>(&dims);
  }
  std::vector<unsigned int> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

typedef Types<float, double, std::complex<float>, std::complex<double>, float_complext, double_complext> Implementations;

TYPED_TEST_CASE(hoNDArray_blas_Test, Implementations);

TYPED_TEST(hoNDArray_blas_Test,dotTest){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(),real(dot(&this->Array,&this->Array)));
  fill(&this->Array2,TypeParam(2));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements()*2,real(dot(&this->Array,&this->Array2)));
}

TYPED_TEST(hoNDArray_blas_Test,axpyTest){
  fill(&this->Array,TypeParam(71));
  fill(&this->Array2,TypeParam(97));
  axpy(TypeParam(11),&this->Array,&this->Array2);
  TypeParam val = this->Array2.get_data_ptr()[10];
  EXPECT_FLOAT_EQ(878,real(val));
}

TYPED_TEST(hoNDArray_blas_Test,nrm2Test){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(std::sqrt((double)this->Array.get_number_of_elements()),nrm2(&this->Array));
  fill(&this->Array,TypeParam(3));
  EXPECT_FLOAT_EQ(std::sqrt(3.0*3.0*this->Array.get_number_of_elements()),nrm2(&this->Array));
}

TYPED_TEST(hoNDArray_blas_Test,asumTest){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements(),real(asum(&this->Array)));
  fill(&this->Array,TypeParam(-3));
  EXPECT_FLOAT_EQ(this->Array.get_number_of_elements()*3,real(asum(&this->Array)));
}

TYPED_TEST(hoNDArray_blas_Test,aminTest){
  fill(&this->Array,TypeParam(100));
  this->Array.get_data_ptr()[23]=TypeParam(-50);
  EXPECT_EQ(23,amin(&this->Array));
  this->Array.get_data_ptr()[48]=TypeParam(2);
  EXPECT_EQ(48,amin(&this->Array));
}

TYPED_TEST(hoNDArray_blas_Test,amaxTest){
  fill(&this->Array,TypeParam(1));
  this->Array.get_data_ptr()[23]=TypeParam(2);
  EXPECT_EQ(23,amax(&this->Array));
  this->Array.get_data_ptr()[48]=TypeParam(-50);
  EXPECT_EQ(48,amax(&this->Array));
}
