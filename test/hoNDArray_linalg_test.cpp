#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"

#include <numeric>
#include <functional>
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

TYPED_TEST_SUITE(hoNDArray_linalg_TestReal, realImplementations);

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

TEST(hoNDArray_linalg, posv_test) {


    std::vector<std::complex<float>> Adata = {{5.96, 0.00},   {0.40, -1.19}, {-0.83, -0.48}, {-0.57, 0.40},
                                              {0.40, 1.19},   {7.95, 0.00},  {0.33, 0.09},   {0.22, 0.74},
                                              {-0.83, 0.48},  {0.33, -0.09}, {4.43, 0.00},   {-1.09, 0.32},
                                              {-0.57, -0.40}, {0.22, -0.74}, {-1.09, -0.32}, {3.46, 0.00}};

    hoNDArray<std::complex<float>> A(4, 4, Adata.data());
    A = permute(A, {1, 0});
   

    //EXPECT_EQ(1,2);
    std::vector<std::complex<float>> Bdata = {{-2.94, 5.79}, {8.44, 3.07},  {8.12, -9.12}, {1.00, -4.62},
                                          {9.09, -5.03}, {3.64, -2.33}, {7.36, 6.77},  {8.04, 2.87}};

    hoNDArray<std::complex<float>> B(2, 4,Bdata.data());

    B = permute(B, {1, 0});


    std::vector<std::complex<float>> solutionData = {{0.80, 1.62},  {2.52, 0.61},  {1.26, -1.78}, {0.01, -1.38},
                                                     {3.38, -0.29}, {2.42, -0.52}, {3.46, 2.92},  {3.77, 1.37}};

    hoNDArray<std::complex<float>> solution(2, 4, solutionData.data());

    solution = permute(solution, {1, 0});

   
    posv(A, B);


    bool all_near = std::inner_product(B.begin(), B.end(), solution.begin(),true, std::logical_and(),
                                       [](auto v1, auto v2) { return abs(v1 - v2) < 0.01; });

    EXPECT_TRUE(all_near);
}


