#include "hoNDArray_linalg.h"


#include <gtest/gtest.h>
#include <complex>
#include <vector>

using namespace Gadgetron;
using testing::Types;

template <typename T> class hoNDArray_linalg_TestReal : public ::testing::Test {
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

TYPED_TEST_CASE(hoNDArray_linalg_TestReal, realImplementations);

TYPED_TEST(hoNDArray_linalg_TestReal, linFitTest)
{
    hoNDArray<TypeParam> x, y;
    x.create(5);

    x(0) = 0;
    x(1) = 1;
    x(2) = 2;
    x(3) = 3;
    x(4) = 4;

    TypeParam a = 2.0, b=4.0;

    y.create(5);

    size_t n;
    for (n=0; n<y.get_number_of_elements(); n++)
    {
        y(n) = a*x(n) + b;
    }

    TypeParam ra, rb;
    Gadgetron::linFit(x, y, ra, rb);

    EXPECT_LT(std::abs(ra - a)/a, 1e-3);
    EXPECT_LT(std::abs(rb - b) / a, 1e-3);
}


