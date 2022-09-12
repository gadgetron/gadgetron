/** \file       curveFitting_test.cpp
    \brief      Test case for curve fitting and its solvers

    \author     Hui Xue
*/

#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"
#include "hoNDHarrWavelet.h"
#include "hoNDRedundantWavelet.h"
#include "hoNDArray_math.h"
#include "simplexLagariaSolver.h"
#include "twoParaExpDecayOperator.h"
#include "twoParaExpRecoveryOperator.h"
#include "curveFittingCostFunction.h"
#include "cmr_t1_mapping.h"
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
TYPED_TEST_SUITE(curveFitting_test, realImplementations);

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

TYPED_TEST(curveFitting_test, T1SR)
{
    typedef Gadgetron::twoParaExpRecoveryOperator< std::vector<TypeParam> > SignalType;
    typedef Gadgetron::leastSquareErrorCostFunction< std::vector<TypeParam> > CostType;

    // define solver
    Gadgetron::simplexLagariaSolver< std::vector<TypeParam>, SignalType, CostType > solver;

    // define signal model
    SignalType t1_sr;

    // define cost function
    CostType lse;

    solver.signal_model_ = &t1_sr;
    solver.cf_ = &lse;

    solver.max_iter_ = 150;
    solver.max_fun_eval_ = 1000;
    solver.thres_fun_ = 1e-4;

    // set measured points
    solver.x_.resize(11, 545); // echo time, in ms
    solver.x_[10] = 10000;

    solver.y_.resize(11); // intensity
    solver.y_[0] = 178;
    solver.y_[1] = 185;
    solver.y_[2] = 182;
    solver.y_[3] = 189;
    solver.y_[4] = 178;
    solver.y_[5] = 180;
    solver.y_[6] = 187;
    solver.y_[7] = 179;
    solver.y_[8] = 177;
    solver.y_[9] = 177;
    solver.y_[10] = 471;

    std::vector<TypeParam> b(2), guess(2, 0);

    guess[0] = *std::max_element(solver.y_.begin(), solver.y_.end());
    guess[1] = solver.x_[solver.x_.size() / 2];

    solver.solve(b, guess);

    EXPECT_NEAR(b[0], 471.062894, 0.003);
    EXPECT_NEAR(b[1], 1122.36963, 0.003);
}

TYPED_TEST(curveFitting_test, T1SRMapping)
{
    Gadgetron::ImageIOAnalyze gt_exporter_;
    Gadgetron::GadgetronTimer gt_timer(false);

    Gadgetron::CmrT1SRMapping<float> t1_sr;

    t1_sr.fill_holes_in_maps_ = true;
    t1_sr.max_size_of_holes_ = 20;
    t1_sr.hole_marking_value_ = 0;
    t1_sr.compute_SD_maps_ = true;

    t1_sr.ti_.resize(11, 545);
    t1_sr.ti_[10] = 10000;

    t1_sr.max_iter_ = 150;
    t1_sr.thres_fun_ = 1e-4;
    t1_sr.max_map_value_ = 4000;

    t1_sr.verbose_ = true;
    // t1_sr.debug_folder_ = debug_folder_full_path_;
    t1_sr.perform_timing_ = true;

    size_t RO = 192;
    size_t E1 = 144;
    size_t N = t1_sr.ti_.size();
    size_t S = 1;
    size_t SLC = 1;

    std::vector<float> y(11);
    y[0] = 178;
    y[1] = 185;
    y[2] = 182;
    y[3] = 189;
    y[4] = 178;
    y[5] = 180;
    y[6] = 187;
    y[7] = 179;
    y[8] = 177;
    y[9] = 177;
    y[10] = 471;

    t1_sr.data_.create(RO, E1, N, S, SLC);

    size_t slc, s, n;

    for (slc = 0; slc < SLC; slc++)
    {
        for (s = 0; s < S; s++)
        {
            Gadgetron::hoNDArray<float> data2D;
            for (n = 0; n < N; n++)
            {
                data2D.create(RO, E1, &(t1_sr.data_(0, 0, n, s, slc)));
                Gadgetron::fill(data2D, y[n]);
            }
        }
    }

    t1_sr.mask_for_mapping_.create(RO, E1, SLC);

    // get the image with longest TS time
    hoNDArray<float> mag_longest_TS;
    mag_longest_TS.create(RO, E1, SLC);

    for (slc = 0; slc < SLC; slc++)
    {
        memcpy(&mag_longest_TS(0, 0, slc), &t1_sr.data_(0, 0, N - 1, 0, slc), sizeof(float)*RO*E1);
    }

    double scale_factor = 1.0;
    GDEBUG_STREAM("CmrParametricT1SRMappingGadget, find incoming image has scale factor of " << scale_factor);

    t1_sr.mask_for_mapping_.create(RO, E1, SLC);
    Gadgetron::fill(t1_sr.mask_for_mapping_, (float)1);

    t1_sr.mask_for_mapping_(12, 23, 0) = 0;
    t1_sr.mask_for_mapping_(12, 24, 0) = 0;
    t1_sr.mask_for_mapping_(13, 23, 0) = 0;

    // -------------------------------------------------------------
    // perform mapping

    { gt_timer.start("CmrParametricT1SRMappingGadget, t1_sr.perform_parametric_mapping"); }
    t1_sr.perform_parametric_mapping();
    { gt_timer.stop(); }

    size_t num_para = t1_sr.get_num_of_paras();

    // -------------------------------------------------------------
    // get the results

    t1_sr.map_.print(std::cout);
    t1_sr.para_.print(std::cout);

    EXPECT_NEAR(t1_sr.para_(0, 0, 0, 0, 0), 471.062894, 0.003);
    EXPECT_NEAR(t1_sr.para_(RO / 2, E1 / 2, 0, 0, 0), 471.062894, 0.003);

    EXPECT_NEAR(t1_sr.map_(0, 0, 0, 0), 1122.36963, 0.003);
    EXPECT_NEAR(t1_sr.map_(RO/2, E1/2, 0, 0), 1122.36963, 0.003);
    EXPECT_NEAR(t1_sr.map_(RO-1, E1-1, 0, 0), 1122.36963, 0.003);
    EXPECT_NEAR(t1_sr.map_(37, 86, 0, 0), 1122.36963, 0.003);

    // test hole filling
    EXPECT_NEAR(t1_sr.map_(12, 23, 0, 0), 1122.36963, 1.0);
}
