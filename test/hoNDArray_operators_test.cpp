#include "hoNDArray_elemwise.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_operators_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    size_t vdims2[] = {37, 49}; //Smaller dimensionality to test batch mode
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    dims2 = std::vector<size_t>(vdims2,vdims2+sizeof(vdims2)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims2);
  }
  std::vector<size_t> dims;
  std::vector<size_t> dims2;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

template <typename T> class hoNDArray_operators_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    size_t vdims2[] = {37, 49}; //Smaller dimensionality to test batch mode
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    dims2 = std::vector<size_t>(vdims2,vdims2+sizeof(vdims2)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
    Array2 = hoNDArray<T>(dims2);
  }
  std::vector<size_t> dims;
  std::vector<size_t> dims2;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

typedef Types<float, double> realImplementations;
typedef Types<std::complex<float>, std::complex<double>, float_complext, double_complext> cplxImplementations;

TYPED_TEST_SUITE(hoNDArray_operators_TestReal, realImplementations);

TYPED_TEST(hoNDArray_operators_TestReal,equalsAddTest1){
  TypeParam v1 = TypeParam(46865.35435);
  TypeParam v2 = TypeParam(13784.34);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array += this->Array2;
  EXPECT_FLOAT_EQ(v1+v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsAddTest2){
  TypeParam v1 = TypeParam(98.4);
  TypeParam v2 = TypeParam(2.2);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array += v2;
  EXPECT_FLOAT_EQ(v1+v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsSubtractTest1){
  TypeParam v1 = TypeParam(98475334.34);
  TypeParam v2 = TypeParam(2452.234);
  unsigned int idx = 124999;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array -= this->Array2;
  EXPECT_FLOAT_EQ(v1-v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsSubtractTest2){
  TypeParam v1 = TypeParam(4.4);
  TypeParam v2 = TypeParam(9212.21);
  unsigned int idx = 122131;
  fill(&this->Array,v1);
  this->Array -= v2;
  EXPECT_FLOAT_EQ(v1-v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsMultiplyTest1){
  TypeParam v1 = TypeParam(342.145);
  TypeParam v2 = TypeParam(43545.43);
  unsigned int idx = 12344;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array *= this->Array2;
  EXPECT_FLOAT_EQ(v1*v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsMultiplyTest2){
  TypeParam v1 = TypeParam(43534.443);
  TypeParam v2 = TypeParam(92.842);
  unsigned int idx = 96735;
  fill(&this->Array,v1);
  this->Array *= v2;
  EXPECT_FLOAT_EQ(v1*v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsDivideTest1){
  TypeParam v1 = TypeParam(644.24);
  TypeParam v2 = TypeParam(38564.64);
  unsigned int idx = 98322;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array /= this->Array2;
  EXPECT_FLOAT_EQ(v1/v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST(hoNDArray_operators_TestReal,equalsDivideTest2){
  TypeParam v1 = TypeParam(56342.24);
  TypeParam v2 = TypeParam(23434.34);
  unsigned int idx = 12591;
  fill(&this->Array,v1);
  this->Array /= v2;
  EXPECT_FLOAT_EQ(v1/v2,this->Array.get_data_ptr()[idx]);
}

TYPED_TEST_SUITE(hoNDArray_operators_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_operators_TestCplx,equalsAddTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array += this->Array2;
  EXPECT_FLOAT_EQ(real(v1+v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1+v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsAddTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array += v2;
  EXPECT_FLOAT_EQ(real(v1+v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1+v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsAddTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array += real(v2);
  EXPECT_FLOAT_EQ(real(v1+v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1+v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsSubtractTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array -= this->Array2;
  EXPECT_FLOAT_EQ(real(v1-v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1-v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsSubtractTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array -= v2;
  EXPECT_FLOAT_EQ(real(v1-v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1-v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsSubtractTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array -= real(v2);
  EXPECT_FLOAT_EQ(real(v1-v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1-v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsMultiplyTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array *= this->Array2;
  EXPECT_FLOAT_EQ(real(v1*v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1*v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsMultiplyTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12925;
  fill(&this->Array,v1);
  this->Array *= v2;
  EXPECT_FLOAT_EQ(real(v1*v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1*v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsMultiplyTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array *= real(v2);
  EXPECT_FLOAT_EQ(real(v1*v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1*v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsDivideTest1){
  TypeParam v1 = TypeParam(46865.35435, 534544.534523);
  TypeParam v2 = TypeParam(13784.34, 54543543.1243);
  unsigned int idx = 73243;
  fill(&this->Array,v1);
  fill(&this->Array2,v2);
  this->Array /= this->Array2;
  EXPECT_FLOAT_EQ(real(v1/v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1/v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsDivideTest2){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,3.23);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array /= v2;
  EXPECT_FLOAT_EQ(real(v1/v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1/v2),imag(this->Array.get_data_ptr()[idx]));
}

TYPED_TEST(hoNDArray_operators_TestCplx,equalsDivideTest3){
  TypeParam v1 = TypeParam(98.4, 45.34);
  TypeParam v2 = TypeParam(2.2,0.0);
  unsigned int idx = 12295;
  fill(&this->Array,v1);
  this->Array /= real(v2);
  EXPECT_FLOAT_EQ(real(v1/v2),real(this->Array.get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(v1/v2),imag(this->Array.get_data_ptr()[idx]));
}

