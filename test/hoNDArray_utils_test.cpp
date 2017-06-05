#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "complext.h"
#include "GadgetronTimer.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_utils_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
};

template <typename T> class hoNDArray_utils_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
};

typedef Types<float, double> realImplementations;
typedef Types</*std::complex<float>, std::complex<double>,*/ float_complext, double_complext> cplxImplementations;

TYPED_TEST_CASE(hoNDArray_utils_TestReal, realImplementations);

TYPED_TEST(hoNDArray_utils_TestReal, fillTest)
{
    hoNDArray<int> src;
    src.create(34, 25, 58, 37);
    Gadgetron::fill(src, 1);

    hoNDArray<int> dst;
    dst.create(55, 60, 78, 87);
    Gadgetron::fill(dst, 2);

    vector_td<size_t, 3> offset_src;
    offset_src[0] = 12;
    offset_src[1] = 4;
    offset_src[2] = 13;

    vector_td<size_t, 3> size;
    size[0] = 8;
    size[1] = 20;
    size[2] = 33;

    vector_td<size_t, 3> offset_dst;
    offset_dst[0] = 0;
    offset_dst[1] = 25;
    offset_dst[2] = 40;

    Gadgetron::fill(offset_src, size, &src, offset_dst, &dst);

    EXPECT_EQ(2, dst(offset_dst[0], offset_dst[1] - 1, offset_dst[2], 0));
    EXPECT_EQ(2, dst(offset_dst[0], offset_dst[1], offset_dst[2] - 1, 0));
    EXPECT_EQ(2, dst(offset_dst[0], offset_dst[1] + size[1], offset_dst[2], 0));
    EXPECT_EQ(2, dst(offset_dst[0], offset_dst[1], offset_dst[2] + size[2], 0));
    EXPECT_EQ(2, dst(offset_dst[0], offset_dst[1] + size[1], offset_dst[2] + size[2], 0));

    EXPECT_EQ(1, dst(offset_dst[0], offset_dst[1], offset_dst[2], 0));
    EXPECT_EQ(1, dst(offset_dst[0], offset_dst[1] + 1, offset_dst[2], 0));
    EXPECT_EQ(1, dst(offset_dst[0], offset_dst[1], offset_dst[2] + 1, 0));
    EXPECT_EQ(1, dst(offset_dst[0], offset_dst[1] + size[1] - 1, offset_dst[2], 0));
    EXPECT_EQ(1, dst(offset_dst[0], offset_dst[1], offset_dst[2] + size[2] - 1, 0));
    EXPECT_EQ(1, dst(offset_dst[0], offset_dst[1] + size[1] - 1, offset_dst[2] + size[2] - 1, 0));
}

TYPED_TEST(hoNDArray_utils_TestReal,permuteTest){

  fill(&this->Array,TypeParam(1));

  std::vector<size_t> order;
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

TYPED_TEST(hoNDArray_utils_TestReal,shiftDimTest){

  fill(&this->Array,TypeParam(1));
  this->Array.get_data_ptr()[37] = 2;

  EXPECT_FLOAT_EQ(1, shift_dim(&this->Array,0)->at(0));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,0)->at(37));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,1)->at(1));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,-1)->at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,2)->at(23*37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,3)->at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(&this->Array,4)->at(37));
}

TYPED_TEST(hoNDArray_utils_TestReal,sumTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(49*v1,sum(&this->Array,1)->get_data_ptr()[idx]);

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(23*v1,sum(&this->Array,2)->get_data_ptr()[idx]);

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(19*v1,sum(&this->Array,3)->get_data_ptr()[idx]);
}

TYPED_TEST_CASE(hoNDArray_utils_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_utils_TestCplx,permuteTest){

  fill(&this->Array,TypeParam(1,1));

  std::vector<size_t> order;
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

  order.clear();
  order.push_back(0); order.push_back(1); order.push_back(3); order.push_back(2);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array, &order)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array, &order)->at(37)));

  order.clear();
  order.push_back(0); order.push_back(2); order.push_back(3); order.push_back(1);

  EXPECT_FLOAT_EQ(2, real(permute(&this->Array, &order)->at(37*23*19)));
  EXPECT_FLOAT_EQ(3, imag(permute(&this->Array, &order)->at(37*23*19)));
}

TYPED_TEST(hoNDArray_utils_TestCplx,shiftDimTest){

  fill(&this->Array,TypeParam(1,1));
  this->Array.get_data_ptr()[37]=TypeParam(2,3);

  EXPECT_FLOAT_EQ(1, real(shift_dim(&this->Array,0)->at(0)));
  EXPECT_FLOAT_EQ(1, imag(shift_dim(&this->Array,0)->at(0)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,0)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,0)->at(37)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,1)->at(1)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,1)->at(1)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,-1)->at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,-1)->at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,2)->at(23*37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,2)->at(23*37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,3)->at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,3)->at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(&this->Array,4)->at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(&this->Array,4)->at(37)));
}

TYPED_TEST(hoNDArray_utils_TestCplx,sumTest){
  TypeParam v1 = TypeParam(12.34, 56.78);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(49)*v1),real(sum(&this->Array,1)->get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(49)*v1),imag(sum(&this->Array,1)->get_data_ptr()[idx]));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(23)*v1),real(sum(&this->Array,2)->get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(23)*v1),imag(sum(&this->Array,2)->get_data_ptr()[idx]));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(19)*v1),real(sum(&this->Array,3)->get_data_ptr()[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(19)*v1),imag(sum(&this->Array,3)->get_data_ptr()[idx]));
}
