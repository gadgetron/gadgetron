#include "gtest/gtest.h"
#include "complext.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoSDC.h"
#include "vector_td_utilities.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

using namespace Gadgetron;
using testing::Types;

template<typename Real>
class hoSDC_test : public ::testing::Test
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
                GDEBUG_STREAM("Unit Test for Gadgetron hoSDC");
                gt_ut_data_folder_ = gt_ut_folder_;
                gt_ut_res_folder_ = gt_ut_folder_;
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

TYPED_TEST_SUITE(hoSDC_test, realImplementations);

TYPED_TEST(hoSDC_test, randomTestOne)
{
    typedef float T;

    if (this->gt_ut_folder_.empty())
    {
        GDEBUG_STREAM("Gadgetron unit test directory is not set ... ");
        return;
    }

    // Load test trajectory.
    hoNDArray<T> k_spiral;
    this->gt_io_.import_array(k_spiral, this->gt_ut_data_folder_ + "/spiral/k_spiral");
    std::stringstream k_spiral_stream;
    k_spiral.print(k_spiral_stream);
    GINFO(k_spiral_stream.str().c_str());
    T v = Gadgetron::nrm2(k_spiral); GDEBUG_STREAM("k_spiral = " << v);

    // Copy test trajectory to array of vector_td elements.
    size_t num_samples = k_spiral.get_number_of_elements() / 2;
    hoNDArray<vector_td<T, 2>> traj;
    traj.create(num_samples);
    for(size_t n = 0; n < num_samples; n++)
    {
        traj(n)[0] = k_spiral(n, 0);
        traj(n)[1] = k_spiral(n, 1);
    }

    // Run density estimation.
    vector_td< size_t, 2 > dims;
    dims[0] = 256;
    dims[1] = 256;

    this->timer_.start("Estimate DCW...");
    hoNDArray<T> res = *estimate_dcw(traj, dims, T(3), 20);
    this->timer_.stop();
    std::stringstream res_stream;
    res.print(res_stream);
    GINFO(res_stream.str().c_str());

    // Save resulting weights.
    this->gt_io_.export_array(res, this->gt_ut_res_folder_ + "/spiral/res_weights");

    // Load reference weights.
    hoNDArray<T> ref;
    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/spiral/ref_weights");
    std::stringstream ref_stream;
    ref.print(ref_stream);
    GINFO(ref_stream.str().c_str());
    v = Gadgetron::nrm2(ref); GDEBUG_STREAM("ref = " << v);

    // Compare resulting weights against reference.
    hoNDArray<T> diff;
    Gadgetron::subtract(res, ref, diff);
    v = Gadgetron::nrm2(diff);
    T norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(v/norm_ref, 0.00001);
}
