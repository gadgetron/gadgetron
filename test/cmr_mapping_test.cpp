#include "gtest/gtest.h"
#include "complext.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNFFT.h"
#include "cmr_t2_mapping.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

using namespace Gadgetron;
using testing::Types;

template<typename Real>
class cmr_mapping_test : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            GDEBUG_STREAM("=============================================================================================");
            char* ut_dir = ::getenv("GT_UNITTEST_DIRECTORY");
            if(ut_dir!=NULL)
            {
                gt_ut_folder_ = std::string(ut_dir);
            }
            else
                gt_ut_folder_.clear();

            GDEBUG_STREAM("=============================================================================================");

            if(!gt_ut_folder_.empty())
            {
                GDEBUG_STREAM("Unit Test for Gadgetron hoNFFT");
                gt_ut_data_folder_ = gt_ut_folder_;
                gt_ut_res_folder_ = gt_ut_folder_ + "/../result/";
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

TYPED_TEST_SUITE(cmr_mapping_test, realImplementations);

TYPED_TEST(cmr_mapping_test, T2Mapping)
{
    typedef float T;

    if (this->gt_ut_folder_.empty())
    {
        GDEBUG_STREAM("Gadgetron unit test directory is not set ... ");
        return;
    }

    hoNDArray< T > data;
    this->gt_io_.import_array(data, this->gt_ut_data_folder_ + "/T2Mapping/MultiEcho_Image_Mag_REP0_SLC0");
    data.print(std::cout);

    T v = Gadgetron::nrm2(data); GDEBUG_STREAM("data = " << v);

    size_t RO = data.get_size(0);
    size_t E1 = data.get_size(1);
    size_t SET = data.get_size(2);

    hoNDArray< T > dataArray;
    dataArray.create(RO, E1, SET, 1, 1, data.get_data_ptr());

    CmrT2Mapping<T> t2mapper;
    t2mapper.verbose_ = true;
    t2mapper.debug_folder_ = this->gt_ut_res_folder_ + "/T2Mapping/";
    t2mapper.perform_timing_ = true;

    t2mapper.compute_SD_maps_ = true;

    t2mapper.data_ = dataArray;
    t2mapper.ti_.resize(3);
    t2mapper.ti_[0] = 0;
    t2mapper.ti_[1] = 25;
    t2mapper.ti_[2] = 55;

    std::vector<T> yi(3);
    yi[0] = 738.33;
    yi[1] = 436.36;
    yi[2] = 244.88;

    std::vector<T> guess, bi, sd;
    t2mapper.get_initial_guess(t2mapper.ti_, yi, guess);

    T t2(0), t2_sd(0);
    t2mapper.compute_map(t2mapper.ti_, yi, guess, bi, t2);

    t2mapper.compute_sd(t2mapper.ti_, yi, bi, sd, t2_sd);

    t2mapper.perform_parametric_mapping();

    this->gt_io_.export_array(t2mapper.map_, this->gt_ut_res_folder_ + "/T2Mapping/t2_map");
    this->gt_io_.export_array(t2mapper.sd_map_, this->gt_ut_res_folder_ + "/T2Mapping/t2_sd_map");

    // compare against ground truth
    hoNDArray<T> ref;
    hoNDArray<T> diff;
    T norm_ref;

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/T2Mapping/ref_t2_map");
    Gadgetron::subtract(t2mapper.map_, ref, diff);
    v = Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(v/norm_ref, 0.002);

    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/T2Mapping/ref_sd_t2_map");
    Gadgetron::subtract(t2mapper.sd_map_, ref, diff);
    v =  Gadgetron::nrm2(diff);
    norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(v / norm_ref, 0.002);
}
