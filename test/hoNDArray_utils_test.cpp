#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "complext.h"
#include "GadgetronTimer.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include <range/v3/view.hpp>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_utils_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
};

template <typename T> class hoNDArray_utils_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(dims);
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

  EXPECT_FLOAT_EQ(1, permute(this->Array,order)[0]);
  EXPECT_FLOAT_EQ(2, permute(this->Array,order)[37]);

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order)[1]);

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order)[19]);

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order)[851]);
}

TYPED_TEST(hoNDArray_utils_TestReal,shiftDimTest){

  fill(&this->Array,TypeParam(1));
  this->Array.get_data_ptr()[37] = 2;

  EXPECT_FLOAT_EQ(1, shift_dim(this->Array,0)[0]);
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,0)[37]);
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,1)[1]);
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,-1)[37*19]);
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,2)[23*37*19]);
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,3)[37*19]);
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,4)[37]);
}

TYPED_TEST(hoNDArray_utils_TestReal,sumTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(49*v1,sum(this->Array,1)[idx]);

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(23*v1,sum(this->Array,2)[idx]);

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(19*v1,sum(this->Array,3)[idx]);
}

TYPED_TEST_CASE(hoNDArray_utils_TestCplx, cplxImplementations);

TYPED_TEST(hoNDArray_utils_TestCplx,permuteTest){

  fill(&this->Array,TypeParam(1,1));

  std::vector<size_t> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);
  
  this->Array.get_data_ptr()[37] = TypeParam(2,3);

  EXPECT_FLOAT_EQ(1, real(permute(this->Array,order)[0]));
  EXPECT_FLOAT_EQ(1, imag(permute(this->Array,order)[0]));

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order)[37]));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order)[37]));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order)[1]));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order)[1]));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order)[19]));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order)[19]));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order)[851]));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order)[851]));

  order.clear();
  order.push_back(0); order.push_back(1); order.push_back(3); order.push_back(2);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array, order)[37]));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array, order)[37]));

  order.clear();
  order.push_back(0); order.push_back(2); order.push_back(3); order.push_back(1);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array, order)[37*23*19]));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array, order)[37*23*19]));
}

TYPED_TEST(hoNDArray_utils_TestCplx,shiftDimTest){

  fill(&this->Array,TypeParam(1,1));
  this->Array.get_data_ptr()[37]=TypeParam(2,3);

  EXPECT_FLOAT_EQ(1, real(shift_dim(this->Array,0)[0]));
  EXPECT_FLOAT_EQ(1, imag(shift_dim(this->Array,0)[0]));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,0)[37]));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,0)[37]));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,1)[1]));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,1)[1]));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,-1)[37*19]));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,-1)[37*19]));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,2)[23*37*19]));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,2)[23*37*19]));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,3)[37*19]));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,3)[37*19]));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,4)[37]));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,4)[37]));
}

TYPED_TEST(hoNDArray_utils_TestCplx,sumTest){
  TypeParam v1 = TypeParam(12.34, 56.78);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(49)*v1),real(sum(this->Array,1)[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(49)*v1),imag(sum(this->Array,1)[idx]));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(23)*v1),real(sum(this->Array,2)[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(23)*v1),imag(sum(this->Array,2)[idx]));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(19)*v1),real(sum(this->Array,3)[idx]));
  EXPECT_FLOAT_EQ(imag(TypeParam(19)*v1),imag(sum(this->Array,3)[idx]));
}

TYPED_TEST(hoNDArray_utils_TestReal,repeatTest){

    fill(&this->Array,TypeParam(2));

    unsigned int repeats = 5;
    auto repeated = repeat(this->Array,repeats);

    auto expected_dims = this->Array.dimensions();
    expected_dims.push_back(repeats);

    EXPECT_EQ(expected_dims,repeated.dimensions());

    EXPECT_FLOAT_EQ(sum(this->Array)*repeats,sum(repeated));


}


TEST(hoNDArray_utils_Test,concat_test){
    hoNDArray<float> arr1(19,7);
    arr1.fill(1);
    auto arr2 = arr1;
    arr2.fill(2);
    auto arr3 = arr1;
    arr3.fill(3);

    using namespace ranges;

    auto concatenated = concat( arr1,arr2,arr3);

    EXPECT_EQ(concatenated(12,3,0),1.0f);
    EXPECT_EQ(concatenated(12,3,1),2.0f);
    EXPECT_EQ(concatenated(12,3,2),3.0f);

}


TEST(hoNDArray_utils_Test,concat_test_vector){
    hoNDArray<float> arr1(19,7);
    arr1.fill(1);
    auto arr2 = arr1;
    arr2.fill(2);
    auto arr3 = arr1;
    arr3.fill(3);


    auto concatenated = concat( std::vector{arr1,arr2,arr3});

    EXPECT_EQ(concatenated(12,3,0),1.0f);
    EXPECT_EQ(concatenated(12,3,1),2.0f);
    EXPECT_EQ(concatenated(12,3,2),3.0f);

}

TEST(hoNDArray_utils_Test,concat_along_test){
    hoNDArray<float> arr1(19,7,3);
    arr1.fill(1);

    hoNDArray<float> arr2(19,7,1);
    arr2.fill(2);

    hoNDArray<float> arr3(19,7,5);
    arr3.fill(3);

    auto concatenated = concat_along_dimension( std::vector{arr1,arr2,arr3},2);

    std::vector<size_t> expected_dims = {19,7,9};

    EXPECT_EQ(expected_dims,concatenated.dimensions());

    EXPECT_EQ(concatenated(12,3,0),1.0f);
    EXPECT_EQ(concatenated(12,3,3),2.0f);
    EXPECT_EQ(concatenated(12,3,6),3.0f);

}

TEST(hoNDArray_utils_Test,concat_along_test_non_contigious){
    hoNDArray<float> arr1(19,3,7);
    arr1.fill(1);

    hoNDArray<float> arr2(19,1,7);
    arr2.fill(2);

    hoNDArray<float> arr3(19,5,7);
    arr3.fill(3);

    auto concatenated = concat_along_dimension( std::vector{arr1,arr2,arr3},1);

    std::vector<size_t> expected_dims = {19,9,7};

    EXPECT_EQ(expected_dims,concatenated.dimensions());

    EXPECT_EQ(concatenated(12,0,3),1.0f);
    EXPECT_EQ(concatenated(1,3,6),2.0f);
    EXPECT_EQ(concatenated(7,6,0),3.0f);

}
