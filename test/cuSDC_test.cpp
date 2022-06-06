#include "gtest/gtest.h"
#include "complext.h"
#include "cuNDArray.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "cuSDC.h"
#include "vector_td_utilities.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

using namespace Gadgetron;
using testing::Types;

template<typename Real>
class cuSDC_test : public ::testing::Test
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

TYPED_TEST_SUITE(cuSDC_test, realImplementations);

TYPED_TEST(cuSDC_test, randomTestOne)
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
    k_spiral.print(std::cout);
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
    vector_td<size_t, 2> dims;
    dims[0] = 256;
    dims[1] = 256;

    cuNDArray<vector_td<T, 2>> d_traj(traj);
    this->timer_.start("Estimate DCW...");
    cuNDArray<T> d_res = *estimate_dcw(d_traj, dims, T(3), 20, float(5.5));
    this->timer_.stop();
    hoNDArray<T> res = *d_res.to_host();
    res.print(std::cout);
    
    // Save resulting weights.
    this->gt_io_.export_array(res, this->gt_ut_res_folder_ + "/spiral/res_weights");

    // Load reference weights.
    hoNDArray<T> ref;
    this->gt_io_.import_array(ref, this->gt_ut_data_folder_ + "/spiral/ref_weights");
    ref.print(std::cout);
    v = Gadgetron::nrm2(ref); GDEBUG_STREAM("ref = " << v);

    // Compare resulting weights against reference.
    hoNDArray<T> diff;
    Gadgetron::subtract(res, ref, diff);
    v = Gadgetron::nrm2(diff);
    T norm_ref = Gadgetron::nrm2(ref);
    EXPECT_LE(v/norm_ref, 0.00001);
}
