#include "NFFTOperator.h"
#include "complext.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_reductions.h"
#include "cuNDArray_utils.h"
#include "cuNFFT.h"

#include "vector_td_utilities.h"
#include <complex>
#include <gtest/gtest.h>
#include <vector>

// Std includes
#include "GadgetronTimer.h"
#include "hoArmadillo.h"
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>

using namespace Gadgetron;
using testing::Types;

template <typename T> class cuNFFT_test : public ::testing::Test {
  protected:
    virtual void SetUp() {
        // Prep for NUFFT_TEST
        RO = 2500;
        INT = 1100;
        CHA = 8;
        xsize_ = 256;
        ysize_ = 256;
        kernel_width_ = 3;
        oversampling_factor_ = 2.1;

        std::vector<size_t> data_dims = {RO, INT, CHA};
        fake_data = cuNDArray<float_complext>(data_dims);
        data_dims.pop_back();
        fake_dcw = cuNDArray<float>(data_dims);

        hoNDArray<vector_td<float, 2>> fake_traj_ho(data_dims);
        vector_td<float, 2> init_val;
        init_val[0] = 0.1f;
        init_val[1] = 0.1f;
        fake_traj_ho.fill(init_val);

        fake_traj = cuNDArray<vector_td<float, 2>>(fake_traj_ho);
        fill(&fake_dcw, 1.0f);
        fill(&fake_data, complext(1.0f, 0.0f));

        unsigned int warp_size = cudaDeviceManager::Instance()->warp_size();

        image_dims_.push_back(xsize_);
        image_dims_.push_back(ysize_);

        image_dims_os_ = uint64d2(
            ((static_cast<size_t>(std::ceil(image_dims_[0] * oversampling_factor_)) + warp_size - 1) / warp_size) *
                warp_size,
            ((static_cast<size_t>(std::ceil(image_dims_[1] * oversampling_factor_)) + warp_size - 1) / warp_size) *
                warp_size); // No oversampling is needed in the z-direction for SOS

        recon_dims = {this->image_dims_[0], this->image_dims_[1], CHA};
    }

    cuNDArray<float_complext> fake_data;
    cuNDArray<float> fake_dcw;
    cuNDArray<vector_td<float, 2>> fake_traj;
    size_t RO, INT, Kz, CHA, xsize_, ysize_;
    float kernel_width_, oversampling_factor_;
    uint64d2 image_dims_os_;
    std::vector<size_t> image_dims_;
    boost::shared_ptr<cuNFFT_plan<float, 2>> nfft_plan_;
    std::vector<size_t> recon_dims;
};

typedef Types<float_complext> cplxImplementations;

TYPED_TEST_SUITE(cuNFFT_test, cplxImplementations);

TYPED_TEST(cuNFFT_test, cuNFFT_ATOMIC) {

    std::vector<size_t> flat_dims = {this->fake_traj.get_number_of_elements()};
    cuNDArray<vector_td<float, 2>> flat_traj(flat_dims, this->fake_traj.get_data_ptr());
    this->nfft_plan_ =
        NFFT<cuNDArray, float, 2>::make_plan(from_std_vector<size_t, 2>(this->image_dims_), this->image_dims_os_,
                                             this->kernel_width_, ConvolutionType::ATOMIC);

    {
        GadgetronTimer timer("Preprocess Atomic");

        this->nfft_plan_->preprocess(flat_traj, NFFT_prep_mode::NC2C);
    }
    auto temp = boost::make_shared<cuNDArray<float_complext>>(this->recon_dims);

    {
        GadgetronTimer timer("Recon Atomic");

        this->nfft_plan_->compute(&this->fake_data, *temp, &this->fake_dcw, NFFT_comp_mode::BACKWARDS_NC2C);
    }
}

TYPED_TEST(cuNFFT_test, cuNFFT_STANDARD) {

    std::vector<size_t> flat_dims = {this->fake_traj.get_number_of_elements()};
    cuNDArray<vector_td<float, 2>> flat_traj(flat_dims, this->fake_traj.get_data_ptr());
    this->nfft_plan_ =
        NFFT<cuNDArray, float, 2>::make_plan(from_std_vector<size_t, 2>(this->image_dims_), this->image_dims_os_,
                                             this->kernel_width_, ConvolutionType::STANDARD);
    {
        GadgetronTimer timer("Preprocess Standard");
        this->nfft_plan_->preprocess(flat_traj, NFFT_prep_mode::NC2C);
    }

    auto temp = boost::make_shared<cuNDArray<float_complext>>(this->recon_dims);
    {
        GadgetronTimer timer("Recon Standard");
        this->nfft_plan_->compute(&this->fake_data, *temp, &this->fake_dcw, NFFT_comp_mode::BACKWARDS_NC2C);
    }
}
