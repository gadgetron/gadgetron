#include "hoNDArray_reductions.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_reductions_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = hoNDArray<T>(&dims);
    Array2 = hoNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  hoNDArray<T> Array;
  hoNDArray<T> Array2;
};

typedef Types<float, double> realImplementations;
TYPED_TEST_CASE(hoNDArray_reductions_TestReal, realImplementations);

TYPED_TEST(hoNDArray_reductions_TestReal, sortTest)
{
    hoNDArray<TypeParam> x, r;
    x.create(5);

    x(0) = 5;
    x(1) = 2;
    x(2) = 3;
    x(3) = 4;
    x(4) = 1;

    std::vector<size_t> ind;
    bool isascending = true;
    Gadgetron::sort(x, r, ind, isascending);

    EXPECT_EQ(r(0) , 1);
    EXPECT_EQ(r(1) , 2);
    EXPECT_EQ(r(2) , 3);
    EXPECT_EQ(r(3) , 4);
    EXPECT_EQ(r(4) , 5);

    EXPECT_EQ(ind[0] , 4);
    EXPECT_EQ(ind[1] , 1);
    EXPECT_EQ(ind[2] , 2);
    EXPECT_EQ(ind[3] , 3);
    EXPECT_EQ(ind[4] , 0);

    isascending = false;
    Gadgetron::sort(x, r, ind, isascending);

    EXPECT_EQ(r(0) , 5);
    EXPECT_EQ(r(1) , 4);
    EXPECT_EQ(r(2) , 3);
    EXPECT_EQ(r(3) , 2);
    EXPECT_EQ(r(4) , 1);

    EXPECT_EQ(ind[0] , 0);
    EXPECT_EQ(ind[1] , 3);
    EXPECT_EQ(ind[2] , 2);
    EXPECT_EQ(ind[3] , 1);
    EXPECT_EQ(ind[4] , 4);
}
