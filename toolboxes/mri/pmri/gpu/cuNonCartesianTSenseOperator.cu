#pragma once

#include "cuNonCartesianTSenseOperator.h"
#include "vector_td_utilities.h"

using namespace Gadgetron;

template<class REAL,unsigned int D>
cuNonCartesianTSenseOperator<REAL,D>::cuNonCartesianTSenseOperator(ConvolutionType conv) : cuSenseOperator<REAL,D>() {

    convolutionType = conv;
    is_preprocessed_ = false;
 }

template <class REAL, unsigned int D>
void cuNonCartesianTSenseOperator<REAL, D>::mult_M(cuNDArray<complext<REAL>>* in, cuNDArray<complext<REAL>>* out,
                                                   bool accumulate) {

    if (!in || !out) {
        throw std::runtime_error("cuNonCartesianSenseOperator::mult_M : 0x0 input/output not accepted");
    }
    if (!in->dimensions_equal(&this->domain_dims_) || !out->dimensions_equal(&this->codomain_dims_)) {
        throw std::runtime_error(
            "cuNonCartesianSenseOperator::mult_H: input/output arrays do not match specified domain/codomains");
    }

    std::vector<size_t> full_dimensions = *this->get_domain_dimensions();
    std::vector<size_t> data_dimensions = *this->get_codomain_dimensions();

    full_dimensions.push_back(this->ncoils_);

    std::iter_swap(full_dimensions.end(), full_dimensions.end() - 1); // swap the coil dimension and time

    full_dimensions.pop_back(); // remove time dimension
    cuNDArray<complext<REAL>> tmp(&full_dimensions);

    std::vector<size_t> slice_dimensions = *this->get_domain_dimensions();
    slice_dimensions.pop_back(); // remove time
    auto stride = std::accumulate(slice_dimensions.begin(), slice_dimensions.end(), 1,
                                  std::multiplies<size_t>()); // product of X,Y,and Z

    for (size_t it = 0; it < shots_per_time_.size(); it++) {
        auto slice_view = cuNDArray<complext<REAL>>(slice_dimensions, in->data() + stride * it);
        this->mult_csm(&slice_view, &tmp);

        data_dimensions.pop_back();                     // remove coil dimension from tmp_data;
        data_dimensions.pop_back();                     // remove interleave
        data_dimensions.push_back(shots_per_time_[it]); // insert correct interleave
        data_dimensions.push_back(this->ncoils_);       // insert coils again

        cuNDArray<complext<REAL>> tmp_data(&data_dimensions);
        //   // Forwards NFFT

        if (accumulate) {
            cuNDArray<complext<REAL>> tmp_out(tmp_data.get_dimensions());
            plan_[it]->compute(&tmp, tmp_out, dcw_[it].get(), NFFT_comp_mode::FORWARDS_C2NC);
            tmp_data += tmp_out;
        }

        else
            plan_[it]->compute(tmp, tmp_data, dcw_[it].get(), NFFT_comp_mode::FORWARDS_C2NC);
        size_t inter_acc = 0;
        if (it > 0)
            inter_acc = std::accumulate(shots_per_time_.begin(), shots_per_time_.begin() + (it - 1), 1,
                                        std::multiplies<size_t>()); // product of X,Y,and Z

        cudaMemcpy(out->get_data_ptr() + inter_acc, tmp_data.get_data_ptr(),
                   tmp_data.get_number_of_elements() * sizeof(complext<REAL>), cudaMemcpyDefault);
    }
}

template <class REAL, unsigned int D>
void cuNonCartesianTSenseOperator<REAL, D>::mult_MH(cuNDArray<complext<REAL>>* in, cuNDArray<complext<REAL>>* out,
                                                    bool accumulate) {

    if (!in || !out) {
        throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH : 0x0 input/output not accepted");
    }

    if (!in->dimensions_equal(&this->codomain_dims_) || !out->dimensions_equal(&this->domain_dims_)) {
        throw std::runtime_error(
            "cuNonCartesianSenseOperator::mult_MH: input/output arrays do not match specified domain/codomains");
    }
    std::vector<size_t> tmp_dimensions = *this->get_domain_dimensions();
    std::vector<size_t> tmp_dimensions_data = *this->get_codomain_dimensions();

    auto RO = in->get_size(0);
    auto E1E2 = in->get_size(1);
    auto CHA = in->get_size(2);

    tmp_dimensions_data.pop_back(); // Remove CH dimension
    tmp_dimensions.pop_back();      // Remove the timeDimension

    tmp_dimensions.push_back(this->ncoils_); // add coil dimension
    cuNDArray<complext<REAL>> tmp(&tmp_dimensions);
    tmp_dimensions.pop_back(); // rm coil dimension

    auto stride_ch = std::accumulate(tmp_dimensions_data.begin(), tmp_dimensions_data.end(), 1,
                                     std::multiplies<size_t>()); // product of X,Y,and Z

    auto stride_out = std::accumulate(tmp_dimensions.begin(), tmp_dimensions.end(), 1,
                                      std::multiplies<size_t>()); // product of X,Y,and Z
    for (size_t it = 0; it < shots_per_time_.size(); it++) {
        size_t inter_acc = 0;
        if (it > 0)
            inter_acc = std::accumulate(shots_per_time_.begin(), shots_per_time_.begin() + (it - 1), 1,
                                        std::multiplies<size_t>()); // product of X,Y,and Z
        for (size_t ich = 0; ich < CHA; ich++) {

            tmp_dimensions_data.pop_back(); // Remove INT dimension
            tmp_dimensions_data.push_back(shots_per_time_[it]);

            auto slice_view = cuNDArray<complext<REAL>>(tmp_dimensions_data, in->data() + stride_ch * ich + inter_acc);

            cuNDArray<complext<REAL>> temp_ch_recon(&tmp_dimensions);

            plan_[it]->compute(&slice_view, temp_ch_recon, dcw_[it].get(), NFFT_comp_mode::BACKWARDS_NC2C);

            cudaMemcpy(tmp.get_data_ptr() + tmp_dimensions[0] * tmp_dimensions[1] * tmp_dimensions[2] * ich,
                       temp_ch_recon.get_data_ptr(),
                       tmp_dimensions[0] * tmp_dimensions[1] * tmp_dimensions[2] * sizeof(complext<REAL>),
                       cudaMemcpyDefault);
        }

        if (!accumulate) {
            clear(out);
        }
        auto slice_view_output = cuNDArray<complext<REAL>>(tmp_dimensions, out->data() + stride_out * it);

        this->mult_csm_conj_sum(&tmp, &slice_view_output);
    }
}

template <class REAL, unsigned int D>
void cuNonCartesianTSenseOperator<REAL, D>::setup(_uint64d matrix_size, _uint64d matrix_size_os, REAL W) {
    for (auto ii = 0; ii < shots_per_time_.size(); ii++)
        plan_.push_back(NFFT<cuNDArray, REAL, D>::make_plan(matrix_size, matrix_size_os, W, convolutionType));
}

template <class REAL, unsigned int D>
void cuNonCartesianTSenseOperator<REAL, D>::preprocess(std::vector<cuNDArray<_reald>>& trajectory) {
    if (&(*trajectory.begin()) == 0x0) {
        throw std::runtime_error("cuNonCartesianSenseOperator: cannot preprocess 0x0 trajectory.");
    }

    boost::shared_ptr<std::vector<size_t>> domain_dims = this->get_domain_dimensions();
    if (domain_dims.get() == 0x0 || domain_dims->size() == 0) {
        throw std::runtime_error("cuNonCartesianSenseOperator::preprocess : operator domain dimensions not set");
    }
    for (auto ii = 0; ii < shots_per_time_.size(); ii++)
        plan_[ii]->preprocess(trajectory[ii], NFFT_prep_mode::ALL);
    is_preprocessed_ = true;
}

template <class REAL, unsigned int D>
void cuNonCartesianTSenseOperator<REAL, D>::set_dcw(std::vector<boost::shared_ptr<cuNDArray<REAL>>> dcw) {
    dcw_ = dcw;
}

template <class REAL, unsigned int D>
void cuNonCartesianTSenseOperator<REAL, D>::set_shots_per_time(std::vector<size_t> shots_per_time) {
    shots_per_time_ = shots_per_time;
}

template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<float, 1>;
template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<float, 2>;
template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<float, 3>;
template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<float, 4>;

template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<double, 1>;
template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<double, 2>;
template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<double, 3>;
template class EXPORTGPUPMRI cuNonCartesianTSenseOperator<double, 4>;
