#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_utils_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
    Array = hoNDArray<T>(&dims);
  }
  std::vector<unsigned int> dims;
  hoNDArray<T> Array;
};

template <typename T> class hoNDArray_utils_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    unsigned int vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<unsigned int>(vdims,vdims+sizeof(vdims)/sizeof(unsigned int));
    Array = hoNDArray<T>(&dims);
  }
  std::vector<unsigned int> dims;
  hoNDArray<T> Array;
};

typedef Types<float, double> realImplementations;
typedef Types</*std::complex<float>, std::complex<double>,*/ float_complext, double_complext> cplxImplementations;

TYPED_TEST_CASE(hoNDArray_utils_TestReal, realImplementations);

TYPED_TEST(hoNDArray_utils_TestReal,permuteTest){

  fill(&this->Array,TypeParam(1));

  std::vector<unsigned int> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);

  this->Array.get_data_ptr()[37] = TypeParam(2);

  EXPECT_FLOAT_EQ(1, permute(&this->Array,&order)->at(0));
  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(37));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(1));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(19));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(&this->Array,&order)->at(851));
}


TYPED_TEST_CASE(hoNDArray_utils_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_utils_TestCplx,permuteTest){

  fill(&this->Array,TypeParam(1,1));

  std::vector<unsigned int> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);
  
  this->Array.get_data_ptr()[37] = TypeParam(2,3);

  EXPECT_FLOAT_EQ(1, real(permute(&this->Array,&order)->at(0)));
  EXPECT_FLOAT_EQ(1, imag(permute(&this->Array,&order)->at(0)));

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(37)));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(1)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(1)));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(19)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(19)));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array,&order)->at(851)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array,&order)->at(851)));
}
