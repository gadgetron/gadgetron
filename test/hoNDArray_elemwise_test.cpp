/*
 * hoNDArray_elemwise_test.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: Dae
 */

#include "hoNDArray_elemwise.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_elemwise_TestReal : public ::testing::Test {
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

template <typename T> class hoNDArray_elemwise_TestCplx : public ::testing::Test {
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

typedef Types<float, double> realImplementations;
typedef Types<std::complex<float>, std::complex<double>, float_complext, double_complext> cplxImplementations;

TYPED_TEST_CASE(hoNDArray_elemwise_TestReal, realImplementations);

TYPED_TEST(hoNDArray_elemwise_TestReal,fillTest){
  fill(&this->Array,TypeParam(1.1));
  EXPECT_FLOAT_EQ(1.1,TypeParam(this->Array.get_data_ptr()[5]));
  fill(&this->Array,TypeParam(27.45));
  EXPECT_FLOAT_EQ(27.45,TypeParam(this->Array.get_data_ptr()[3242]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clearTest){
  fill(&this->Array,TypeParam(1));
  EXPECT_FLOAT_EQ(1,TypeParam(this->Array.get_data_ptr()[5324]));
  clear(&this->Array);
  EXPECT_FLOAT_EQ(0,TypeParam(this->Array.get_data_ptr()[5324]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,absTest){
  fill(&this->Array,TypeParam(-5.5));
  EXPECT_FLOAT_EQ(TypeParam(-5.5),TypeParam(this->Array.get_data_ptr()[13]));
  EXPECT_FLOAT_EQ(TypeParam(5.5),TypeParam(abs(&this->Array)->get_data_ptr()[13]));
  fill(&this->Array,TypeParam(-1.3));
  EXPECT_FLOAT_EQ(TypeParam(-1.3),TypeParam(this->Array.get_data_ptr()[2454]));
  abs_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1.3),TypeParam(this->Array.get_data_ptr()[2454]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,sqrtTest){
  fill(&this->Array,TypeParam(17.9));
  EXPECT_FLOAT_EQ(std::sqrt(TypeParam(17.9)),TypeParam(sqrt(&this->Array)->get_data_ptr()[23433]));
  fill(&this->Array,TypeParam(3.14));
  sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(std::sqrt(TypeParam(3.14)),TypeParam(this->Array.get_data_ptr()[32343]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,squareTest){
  fill(&this->Array,TypeParam(1.7));
  EXPECT_FLOAT_EQ(TypeParam(1.7)*TypeParam(1.7),TypeParam(square(&this->Array)->get_data_ptr()[22542]));
  fill(&this->Array,TypeParam(31.4));
  square_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(31.4)*TypeParam(31.4),TypeParam(this->Array.get_data_ptr()[652252]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,reciprocalTest){
  fill(&this->Array,TypeParam(11.7));
  EXPECT_FLOAT_EQ(TypeParam(1)/TypeParam(11.7),TypeParam(reciprocal(&this->Array)->get_data_ptr()[45452]));
  fill(&this->Array,TypeParam(314.114));
  reciprocal_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1)/TypeParam(314.114),TypeParam(this->Array.get_data_ptr()[43432]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,reciprocal_sqrtTest){
  fill(&this->Array,TypeParam(1.9));
  EXPECT_FLOAT_EQ(TypeParam(1)/std::sqrt(TypeParam(1.9)),TypeParam(reciprocal_sqrt(&this->Array)->get_data_ptr()[12345]));
  fill(&this->Array,TypeParam(1.14));
  reciprocal_sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(1)/std::sqrt(TypeParam(1.14)),TypeParam(this->Array.get_data_ptr()[0]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,sgnTest){
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[91] = TypeParam(101.1);
  this->Array.get_data_ptr()[19100] = TypeParam(0);
  EXPECT_FLOAT_EQ(TypeParam(-1),TypeParam(sgn(&this->Array)->get_data_ptr()[28]));
  EXPECT_FLOAT_EQ(TypeParam(1),TypeParam(sgn(&this->Array)->get_data_ptr()[91]));
  EXPECT_FLOAT_EQ(TypeParam(0),TypeParam(sgn(&this->Array)->get_data_ptr()[19100]));
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[9100] = TypeParam(101.1);
  this->Array.get_data_ptr()[19100] = TypeParam(0);
  sgn_inplace(&this->Array);
  EXPECT_FLOAT_EQ(TypeParam(-1),TypeParam(this->Array.get_data_ptr()[2800]));
  EXPECT_FLOAT_EQ(TypeParam(1),TypeParam(this->Array.get_data_ptr()[9100]));
  EXPECT_FLOAT_EQ(TypeParam(0),TypeParam(this->Array.get_data_ptr()[19100]));
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clampTest){
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[354222] = TypeParam(101.3);
  clamp(&this->Array,TypeParam(4.9),TypeParam(100.0));
  EXPECT_FLOAT_EQ(TypeParam(4.9),this->Array.get_data_ptr()[3435]);
  EXPECT_FLOAT_EQ(TypeParam(100.0),this->Array.get_data_ptr()[354222]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clamp_minTest){
  fill(&this->Array,TypeParam(-5.7));
  this->Array.get_data_ptr()[91] = TypeParam(-101.3);
  clamp_min(&this->Array,TypeParam(-10.6));
  EXPECT_FLOAT_EQ(TypeParam(-5.7),this->Array.get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(TypeParam(-10.6),this->Array.get_data_ptr()[91]);
}

TYPED_TEST(hoNDArray_elemwise_TestReal,clamp_maxTest){
  fill(&this->Array,TypeParam(5.7));
  this->Array.get_data_ptr()[91] = TypeParam(101.3);
  clamp_max(&this->Array,TypeParam(10.6));
  EXPECT_FLOAT_EQ(TypeParam(5.7),this->Array.get_data_ptr()[28]);
  EXPECT_FLOAT_EQ(TypeParam(10.6),this->Array.get_data_ptr()[91]);
}

TYPED_TEST_CASE(hoNDArray_elemwise_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_elemwise_TestCplx,fillTest){
  fill(&this->Array,TypeParam(1.1,2.2));
  EXPECT_FLOAT_EQ(1.1,real(TypeParam(this->Array.get_data_ptr()[52323])));
  EXPECT_FLOAT_EQ(2.2,imag(TypeParam(this->Array.get_data_ptr()[52323])));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,clearTest){
  fill(&this->Array,TypeParam(1,1));
  clear(&this->Array);
  EXPECT_FLOAT_EQ(0,real(TypeParam(this->Array.get_data_ptr()[325])));
  EXPECT_FLOAT_EQ(0,imag(TypeParam(this->Array.get_data_ptr()[325])));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,absTest){
  fill(&this->Array,TypeParam(-5.5,7.7));
  EXPECT_FLOAT_EQ(std::sqrt(5.5*5.5+7.7*7.7),abs(&this->Array)->get_data_ptr()[32113]);
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,sqrtTest){
  fill(&this->Array,TypeParam(17.9,3.5));
  EXPECT_FLOAT_EQ(real(sqrt(TypeParam(17.9,3.5))),real(sqrt(&this->Array)->get_data_ptr()[2131]));
  EXPECT_FLOAT_EQ(imag(sqrt(TypeParam(17.9,3.5))),imag(sqrt(&this->Array)->get_data_ptr()[2131]));
  fill(&this->Array,TypeParam(3.14,4.13));
  sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(sqrt(TypeParam(3.14,4.13))),real(this->Array.get_data_ptr()[120000]));
  EXPECT_FLOAT_EQ(imag(sqrt(TypeParam(3.14,4.13))),imag(this->Array.get_data_ptr()[120000]));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,squareTest){
  fill(&this->Array,TypeParam(1.7,7.1));
  EXPECT_FLOAT_EQ(real(TypeParam(1.7,7.1)*TypeParam(1.7,7.1)),real(square(&this->Array)->get_data_ptr()[22123]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1.7,7.1)*TypeParam(1.7,7.1)),imag(square(&this->Array)->get_data_ptr()[22123]));
  fill(&this->Array,TypeParam(31.4,4.31));
  square_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(31.4,4.31)*TypeParam(31.4,4.31)),real(this->Array.get_data_ptr()[51234]));
  EXPECT_FLOAT_EQ(imag(TypeParam(31.4,4.31)*TypeParam(31.4,4.31)),imag(this->Array.get_data_ptr()[51234]));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,reciprocalTest){
  fill(&this->Array,TypeParam(1.9,2.7));
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/TypeParam(1.9,2.7)),real(reciprocal(&this->Array)->get_data_ptr()[11232]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/TypeParam(1.9,2.7)),imag(reciprocal(&this->Array)->get_data_ptr()[11232]));
  fill(&this->Array,TypeParam(1.14,4.32));
  reciprocal_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/TypeParam(1.14,4.32)),real(this->Array.get_data_ptr()[10]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/TypeParam(1.14,4.32)),imag(this->Array.get_data_ptr()[10]));
}

TYPED_TEST(hoNDArray_elemwise_TestCplx,reciprocal_sqrtTest){
  fill(&this->Array,TypeParam(1.9,2.7));
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/sqrt(TypeParam(1.9,2.7))),real(reciprocal_sqrt(&this->Array)->get_data_ptr()[12543]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/sqrt(TypeParam(1.9,2.7))),imag(reciprocal_sqrt(&this->Array)->get_data_ptr()[12543]));
  fill(&this->Array,TypeParam(1.14,4.32));
  reciprocal_sqrt_inplace(&this->Array);
  EXPECT_FLOAT_EQ(real(TypeParam(1,0)/sqrt(TypeParam(1.14,4.32))),real(this->Array.get_data_ptr()[10000]));
  EXPECT_FLOAT_EQ(imag(TypeParam(1,0)/sqrt(TypeParam(1.14,4.32))),imag(this->Array.get_data_ptr()[10000]));
}
