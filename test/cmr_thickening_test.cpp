#include "gtest/gtest.h"
#include "complext.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNFFT.h"
#include "cmr_radial_thickening.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

using namespace Gadgetron;
using testing::Types;

template<typename Real>
class cmr_thickening_test : public ::testing::Test
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
            GDEBUG_STREAM("Unit Test for Gadgetron Radial Thickening");
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

TYPED_TEST_SUITE(cmr_thickening_test, realImplementations);

TYPED_TEST(cmr_thickening_test, 229201050_229201055_362_20190718_retro_cine)
{
    typedef float T;

    if (this->gt_ut_folder_.empty())
    {
        GDEBUG_STREAM("Gadgetron unit test directory is not set ... ");
        return;
    }

    hoNDArray< T > epi_mask;
    this->gt_io_.import_array(epi_mask, this->gt_ut_data_folder_ + "/RadialThickening/retro_cine_epi_mask");
    std::stringstream epi_mask_stream;
    epi_mask.print(epi_mask_stream);
    GINFO(epi_mask_stream.str().c_str());

    hoNDArray< T > endo_mask;
    this->gt_io_.import_array(endo_mask, this->gt_ut_data_folder_ + "/RadialThickening/retro_cine_endo_mask");
    std::stringstream endo_mask_stream;
    endo_mask.print(endo_mask_stream);
    GINFO(endo_mask_stream.str().c_str());

    hoNDArray<T> rad_strain, edge_endo, edge_epi;
    size_t ref_phase = 1;
    Gadgetron::compute_thickening(endo_mask, epi_mask, ref_phase, edge_endo, edge_epi, rad_strain);

    this->gt_io_.export_array(rad_strain, this->gt_ut_res_folder_ + "/RadialThickening/retro_cine_strain");
    this->gt_io_.export_array(edge_endo, this->gt_ut_res_folder_ + "/RadialThickening/retro_edge_endo");
    this->gt_io_.export_array(edge_epi, this->gt_ut_res_folder_ + "/RadialThickening/retro_edge_epi");

    // compare against ground truth
    hoNDArray<T> ref;
    hoNDArray<T> diff;
    T norm_ref;

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/RadialThickening/retro_cine_strain");
    GINFO_STREAM(Gadgetron::nrm2(rad_strain) << std::endl);
    GINFO_STREAM(Gadgetron::nrm2(ref) << std::endl);
    Gadgetron::subtract(ref, rad_strain, diff);
    T v = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(v / norm_ref, 0.1);

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/RadialThickening/retro_cine_edge_endo");
    GINFO_STREAM(Gadgetron::nrm2(edge_endo) << std::endl);
    GINFO_STREAM(Gadgetron::nrm2(ref) << std::endl);
    Gadgetron::subtract(ref, edge_endo, diff);
    T q = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(q / norm_ref, 0.15);

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/RadialThickening/retro_cine_edge_epi");
    GINFO_STREAM(Gadgetron::nrm2(edge_epi) << std::endl);
    GINFO_STREAM(Gadgetron::nrm2(ref) << std::endl);
    Gadgetron::subtract(ref, edge_epi, diff);
    T s = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(s / norm_ref, 0.15);
}
