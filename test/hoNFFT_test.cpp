#include "gtest/gtest.h"
#include "complext.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNFFT.h"
#include "vector_td_utilities.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

using namespace Gadgetron;
using testing::Types;

template<typename Real>
class hoNFFT_2D_NC2C_BACKWARDS : public ::testing::Test{
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

TYPED_TEST_SUITE(hoNFFT_2D_NC2C_BACKWARDS, realImplementations);

TYPED_TEST(hoNFFT_2D_NC2C_BACKWARDS, randomTestOne)
{
    typedef float T;

    if (this->gt_ut_folder_.empty())
    {
        GDEBUG_STREAM("Gadgetron unit test directory is not set ... ");
        return;
    }

    hoNDArray< std::complex<T> > data;
    this->gt_io_.import_array_complex(data, this->gt_ut_data_folder_ + "/spiral/data_spiral_real",
        this->gt_ut_data_folder_ + "/spiral/data_spiral_imag");
    std::stringstream data_stream;
    data.print(data_stream);
    GINFO(data_stream.str().c_str());


    T v = Gadgetron::nrm2(data); GDEBUG_STREAM("data = " << v);

    hoNDArray< T > k_spiral;
    this->gt_io_.import_array(k_spiral, this->gt_ut_data_folder_ + "/spiral/k_spiral");
    std::stringstream k_spiral_stream;
    k_spiral.print(k_spiral_stream);
    GINFO(k_spiral_stream.str().c_str());

    v = Gadgetron::nrm2(k_spiral); GDEBUG_STREAM("k_spiral = " << v);

    size_t num = data.get_size(0);
    size_t CHA = data.get_size(1);

    hoNDArray<vector_td<T, 2>> traj;
    traj.create(num);

    size_t n;
    for(n=0; n<num; n++)
    {
        traj(n)[0] = k_spiral(n, 0);
        traj(n)[1] = k_spiral(n, 1);
    }

    hoNDArray< T > w_spiral;
    this->gt_io_.import_array(w_spiral, this->gt_ut_data_folder_ + "/spiral/w_spiral");
    std::stringstream w_spiral_stream;
    w_spiral.print(w_spiral_stream);
    GINFO(w_spiral_stream.str().c_str());

    v = Gadgetron::nrm2(w_spiral); GDEBUG_STREAM("w_spiral = " << v);

    // test the simple gridding

    hoNDArray< std::complex<T> > res;

    vector_td< size_t, 2 > dims;
    dims[0] = 256;
    dims[1] = 256;

    res.create(2* dims[0], 2* dims[1], CHA);
    Gadgetron::clear(res);

    hoNFFT_plan<T, 2> plan(dims, dims*size_t(2), 3.0);
    plan.preprocess(traj);

    hoNDArray< std::complex<T> > dataCha;
    hoNDArray< std::complex<T> > resCha;

    this->timer_.start("Simple regridding ... ");
    for (size_t cha=0; cha<CHA; cha++)
    {
        dataCha.create(num, data.begin()+cha*num);
        resCha.create(res.get_size(0), res.get_size(1), res.begin()+cha*res.get_size(0)*res.get_size(1));

        plan.compute(dataCha, resCha, &w_spiral,NFFT_comp_mode::BACKWARDS_NC2C);
    }
    this->timer_.stop();

    std::stringstream res_stream;
    res.print(res_stream);
    GINFO(res_stream.str().c_str());


    this->gt_io_.export_array_complex(res, this->gt_ut_res_folder_ + "/spiral/res_regridding");

    // compare against ground truth
    hoNDArray< std::complex<T> > ref;
    this->gt_io_.import_array_complex(ref, this->gt_ut_data_folder_ + "/spiral/ref_regridding_REAL",
        this->gt_ut_data_folder_ + "/spiral/ref_regridding_IMAG");
    std::stringstream ref_stream;
    ref.print(ref_stream);
    GINFO(ref_stream.str().c_str());


    hoNDArray< std::complex<T> > diff;
    Gadgetron::subtract(res, ref, diff);

    v = Gadgetron::nrm2(diff);

    T norm_ref = Gadgetron::nrm2(ref);

    EXPECT_LE(v/norm_ref, 0.00001);
}
