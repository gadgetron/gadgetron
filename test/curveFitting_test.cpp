/** \file       curveFitting_test.cpp
    \brief      Test case for curve fitting and its solvers

    \author     Hui Xue
*/

#include "hoNDHarrWavelet.h"
#include "hoNDRedundantWavelet.h"
#include "hoNDArray_math.h"
#include "simplexLagariaSolver.h"
#include "twoParaExpDecayOperator.h"
#include "curveFittingCostFunction.h"
#include <gtest/gtest.h>
#include <boost/random.hpp>

using namespace Gadgetron;
using testing::Types;

template<typename T> class curveFitting_test : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
    }
};

typedef Types<float> realImplementations;
TYPED_TEST_CASE(curveFitting_test, realImplementations);

TYPED_TEST(curveFitting_test, T2SE)
{
    // define solver
    Gadgetron::simplexLagariaSolver< std::vector<TypeParam>, Gadgetron::twoParaExpDecayOperator< std::vector<TypeParam> >, Gadgetron::leastSquareErrorCostFunction< std::vector<TypeParam> > > solver;

    // define signal model
    Gadgetron::twoParaExpDecayOperator< std::vector<TypeParam> > t2;

    // define cost function
    Gadgetron::leastSquareErrorCostFunction< std::vector<TypeParam> > lse;

    // set measured points
    solver.x_.resize(8); // echo time, in ms
    solver.x_[0] = 10;
    solver.x_[1] = 20;
    solver.x_[2] = 30;
    solver.x_[3] = 40;
    solver.x_[4] = 60;
    solver.x_[5] = 80;
    solver.x_[6] = 120;
    solver.x_[7] = 160;

    solver.y_.resize(8); // intensity
    solver.y_[0] = 606.248226950355;
    solver.y_[1] = 598.40425531914;
    solver.y_[2] = 589.368794326241;
    solver.y_[3] = 580.815602836879;
    solver.y_[4] = 563.170212765957;
    solver.y_[5] = 545.893617021277;
    solver.y_[6] = 512.31914893617;
    solver.y_[7] = 480.723404255319;

    std::vector<TypeParam> bi(2), b(2);
    bi[0] = solver.y_[0];
    bi[1] = 640;

    solver.signal_model_ = &t2;
    solver.cf_ = &lse;

    solver.solve(b, bi);

    EXPECT_NEAR(b[0], 617.257, 0.001);
    EXPECT_NEAR(b[1], 644.417, 0.001);
}
