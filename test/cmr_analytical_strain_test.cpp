#include "gtest/gtest.h"
#include "complext.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNFFT.h"
#include "cmr_analytical_strain.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

using namespace Gadgetron;
using testing::Types;

template<typename Real>
class cmr_analytical_strain_test : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        GDEBUG_STREAM("=============================================================================================");
        char* ut_dir = ::getenv("GT_UNITTEST_DIRECTORY");
        if (ut_dir != NULL)
        {
            gt_ut_folder_ = std::string(ut_dir);
        }
        else
            gt_ut_folder_.clear();

        GDEBUG_STREAM("=============================================================================================");

        if (!gt_ut_folder_.empty())
        {
            GDEBUG_STREAM("Unit Test for Gadgetron Analytical Strain");
            gt_ut_data_folder_ = gt_ut_folder_;
            gt_ut_res_folder_ = gt_ut_folder_ + "/result/";
            GDEBUG_STREAM("gt_ut_data_folder_ is " << gt_ut_data_folder_);
            GDEBUG_STREAM("gt_ut_res_folder_ is " << gt_ut_res_folder_);
        }

        this->timer_.set_timing_in_destruction(false);
    }

    std::string gt_ut_folder_;
    std::string gt_ut_data_folder_;
    std::string gt_ut_res_folder_;

    ImageIOAnalyze gt_io_;
    GadgetronTimer timer_;
};

typedef Types<float> realImplementations;

TYPED_TEST_SUITE(cmr_analytical_strain_test, realImplementations);

TYPED_TEST(cmr_analytical_strain_test, QuadarticStrain)
{
    typedef float T;

    if (this->gt_ut_folder_.empty())
    {
        GDEBUG_STREAM("Gadgetron unit test directory is not set ... ");
        return;
    }

    hoNDArray< T > dx;
    this->gt_io_.import_array(dx, this->gt_ut_data_folder_ + "/SynthesisStrain/grtr_dx");
    dx.print(std::cout);

    hoNDArray< T > dy;
    this->gt_io_.import_array(dy, this->gt_ut_data_folder_ + "/SynthesisStrain/grtr_dy");
    dy.print(std::cout);

    hoNDArray< T > mask;
    this->gt_io_.import_array(mask, this->gt_ut_data_folder_ + "/SynthesisStrain/mask");
    T norm_mask = Gadgetron::nrm2(mask);
    GDEBUG_STREAM("mask is " << norm_mask);

    hoNDArray<double> dx_used, dy_used;

    dx_used.copyFrom(dx);
    dy_used.copyFrom(dy);

    hoNDArray<T> radial, circ;
    this->timer_.start("Compute analytical straion, quad ... ");
    Gadgetron::compute_analytical_strain(dx_used, dy_used, mask, radial, circ);
    this->timer_.stop();

    this->gt_io_.export_array(radial, this->gt_ut_res_folder_ + "/AnalyticalStrain/analytical_quad_rad_strain");
    this->gt_io_.export_array(circ, this->gt_ut_res_folder_ + "/AnalyticalStrain/analytical_quad_circ_strain");

    // compare against ground truth
    hoNDArray<T> ref;
    hoNDArray<T> diff;
    T norm_ref;

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/AnalyticalStrain/analytical_quad_rad_strain");
    std::cout << Gadgetron::nrm2(radial) << std::endl;
    std::cout << Gadgetron::nrm2(ref) << std::endl;
    Gadgetron::subtract(ref, radial, diff);
    T v = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(v / norm_ref, 0.002);

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/AnalyticalStrain/analytical_quad_circ_strain");
    std::cout << Gadgetron::nrm2(circ) << std::endl;
    std::cout << Gadgetron::nrm2(ref) << std::endl;
    Gadgetron::subtract(ref, circ, diff);
    T q = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(q / norm_ref, 0.002);
}

TYPED_TEST(cmr_analytical_strain_test, ConstStrain)
{
    typedef float T;

    if (this->gt_ut_folder_.empty())
    {
        GDEBUG_STREAM("Gadgetron unit test directory is not set ... ");
        return;
    }

    hoNDArray< T > dx;
    this->gt_io_.import_array(dx, this->gt_ut_data_folder_ + "/SynthesisStrain/cons_grtr_dx");
    dx.print(std::cout);

    hoNDArray< T > dy;
    this->gt_io_.import_array(dy, this->gt_ut_data_folder_ + "/SynthesisStrain/cons_grtr_dy");
    dy.print(std::cout);

    hoNDArray< T > mask;
    this->gt_io_.import_array(mask, this->gt_ut_data_folder_ + "/SynthesisStrain/cons_mask");
    T norm_mask = Gadgetron::nrm2(mask);
    GDEBUG_STREAM("mask is " << norm_mask);

    hoNDArray<double> dx_used, dy_used;

    dx_used.copyFrom(dx);
    dy_used.copyFrom(dy);

    hoNDArray<T> radial, circ;
    this->timer_.start("Compute analytical straion, const ... ");
    Gadgetron::compute_analytical_strain(dx_used, dy_used, mask, radial, circ);
    this->timer_.stop();

    this->gt_io_.export_array(radial, this->gt_ut_res_folder_ + "/AnalyticalStrain/analytical_cons_rad_strain");
    this->gt_io_.export_array(circ, this->gt_ut_res_folder_ + "/AnalyticalStrain/analytical_cons_circ_strain");

    // compare against ground truth
    hoNDArray<T> ref;
    hoNDArray<T> diff;
    T norm_ref;

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/AnalyticalStrain/analytical_cons_rad_strain");
    std::cout << Gadgetron::nrm2(radial) << std::endl;
    std::cout << Gadgetron::nrm2(ref) << std::endl;
    Gadgetron::subtract(ref, radial, diff);
    T v = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(v / norm_ref, 0.002);

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/AnalyticalStrain/analytical_cons_circ_strain");
    std::cout << Gadgetron::nrm2(circ) << std::endl;
    std::cout << Gadgetron::nrm2(ref) << std::endl;
    Gadgetron::subtract(ref, circ, diff);
    T q = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(q / norm_ref, 0.002);
}
