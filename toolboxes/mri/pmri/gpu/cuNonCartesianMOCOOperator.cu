#pragma once

#include "cuNonCartesianMOCOOperator.h"
#include "vector_td_utilities.h"

using namespace Gadgetron;

template <class REAL, unsigned int D>
cuNonCartesianMOCOOperator<REAL, D>::cuNonCartesianMOCOOperator(ConvolutionType conv) : cuSenseOperator<REAL, D>() {

    convolutionType = conv;
    is_preprocessed_ = false;
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::mult_M(cuNDArray<complext<REAL>>* in, cuNDArray<complext<REAL>>* out,
                                                 bool accumulate) {

    if (!in || !out) {
        throw std::runtime_error("cuNonCartesianSenseOperator::mult_M : 0x0 input/output not accepted");
    }
    if (!in->dimensions_equal(&this->domain_dims_) || !out->dimensions_equal(&this->codomain_dims_)) {
        throw std::runtime_error(
            "cuNonCartesianSenseOperator::mult_H: input/output arrays do not match specified domain/codomains");
    }
    // Cart -> noncart
    std::vector<size_t> full_dimensions = *this->get_domain_dimensions();   // cart
    std::vector<size_t> data_dimensions = *this->get_codomain_dimensions(); // Non-cart
    data_dimensions.pop_back();                                             // remove coil dimension from tmp_data;

    auto timeD = full_dimensions[full_dimensions.size() - 1];
    full_dimensions.pop_back();
    full_dimensions.push_back(this->ncoils_);
    full_dimensions.push_back(timeD);

    // std::iter_swap(full_dimensions.end(), full_dimensions.end() - 1); // swap the coil dimension and time

    full_dimensions.pop_back(); // remove time dimension

    std::vector<size_t> slice_dimensions = *this->get_domain_dimensions();
    slice_dimensions.pop_back(); // remove time
    auto stride = std::accumulate(slice_dimensions.begin(), slice_dimensions.end(), 1,
                                  std::multiplies<size_t>()); // product of X,Y,and Z

    std::vector<size_t> tmp_dims = *this->get_codomain_dimensions();
    auto stride_data = std::accumulate(tmp_dims.begin(), tmp_dims.end() - 1, 1, std::multiplies<size_t>());
    
    auto input = cuNDArray<complext<REAL>>(slice_dimensions, in->data());

    slice_dimensions.push_back(shots_per_time_.size()); // add time dimenstion
    cuNDArray<complext<REAL>> moving_images(slice_dimensions);
    slice_dimensions.pop_back();
    for (size_t it = 0; it < shots_per_time_.size(); it++) {

        auto inter_acc = std::accumulate(shots_per_time_.begin(), shots_per_time_.begin() + it, size_t(0)) *
                         tmp_dims[0]; // sum of cum sum shots per time


        auto slice_view_in = cuNDArray<complext<REAL>>(slice_dimensions, moving_images.data() + stride * it);
        
        // Move the image to moving image
        slice_view_in = input;
        applyDeformation(&slice_view_in,backward_deformation_[it]);

        cuNDArray<complext<REAL>> tmp(&full_dimensions);
        this->mult_csm(&slice_view_in, &tmp);

        data_dimensions.pop_back();                     // remove interleave
        data_dimensions.push_back(shots_per_time_[it]); // insert correct interleave
        // data_dimensions.push_back(this->ncoils_);       // insert coils again

        // cuNDArray<complext<REAL>> tmp_data(&data_dimensions);

        full_dimensions.pop_back(); // remove ch
        for (size_t iCHA = 0; iCHA < this->ncoils_; iCHA++) {
            auto tmp_view = cuNDArray<complext<REAL>>(full_dimensions, tmp.data() + stride * iCHA);
            auto slice_view_out =
                cuNDArray<complext<REAL>>(data_dimensions, out->data() + inter_acc + stride_data * iCHA);

            if (accumulate) {
                cuNDArray<complext<REAL>> tmp_out(&full_dimensions);
                plan_[it]->compute(tmp_view, tmp_out, &dcw_[it], NFFT_comp_mode::FORWARDS_C2NC);
                slice_view_out += tmp_out;
            } else
                plan_[it]->compute(tmp_view, slice_view_out, &dcw_[it], NFFT_comp_mode::FORWARDS_C2NC);
        }
        full_dimensions.push_back(this->ncoils_);
        // size_t inter_acc = 0;
        // if (it > 0)

        // This is not correct yet ! -- AJ
        // for (size_t iCHA = 0; iCHA < this->ncoils_; iCHA++)
        //     cudaMemcpy(out->get_data_ptr() + inter_acc + stride_data * iCHA,
        //                tmp_data.get_data_ptr() + tmp_data.get_size(0) * tmp_data.get_size(1) * iCHA,
        //                tmp_data.get_size(0) * tmp_data.get_size(1) * sizeof(complext<REAL>), cudaMemcpyDefault);
    }
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::mult_MH(cuNDArray<complext<REAL>>* in, cuNDArray<complext<REAL>>* out,
                                                  bool accumulate) {

    if (!in || !out) {
        throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH : 0x0 input/output not accepted");
    }

    if (!in->dimensions_equal(&this->codomain_dims_) || !out->dimensions_equal(&this->domain_dims_)) {
        throw std::runtime_error(
            "cuNonCartesianSenseOperator::mult_MH: input/output arrays do not match specified domain/codomains");
    }
    std::vector<size_t> out_dimensions = *this->get_domain_dimensions();
    std::vector<size_t> in_dimensions = *this->get_codomain_dimensions();

    auto RO = in->get_size(0);
    auto E1E2 = in->get_size(1);
    auto CHA = in->get_size(2);

    in_dimensions.pop_back(); // Remove CH dimension

    out_dimensions.pop_back();               // Remove the timeDimension
    out_dimensions.push_back(this->ncoils_); // add coil dimension
    cuNDArray<complext<REAL>> tmp(&out_dimensions);
    out_dimensions.pop_back(); // rm coil dimension
    // cuNDArray<complext<REAL>> tmp_coilCmb(&out_dimensions);

    auto stride_ch = std::accumulate(in_dimensions.begin(), in_dimensions.end(), 1,
                                     std::multiplies<size_t>()); // product of X,Y,and Z

    auto stride_out = std::accumulate(out_dimensions.begin(), out_dimensions.end(), 1,
                                      std::multiplies<size_t>()); // product of X,Y,and Z
    if (!accumulate) {
        clear(out);
    }
    out_dimensions.push_back(shots_per_time_.size()); // add time dimension
    cuNDArray<complext<REAL>> moving_images(out_dimensions);
    out_dimensions.pop_back(); // rm time

    auto output = cuNDArray<complext<REAL>>(out_dimensions, out->data() );
    fill<complext<REAL>>(&output,complext<REAL>((REAL)0,(REAL)0));
    for (size_t it = 0; it < shots_per_time_.size(); it++) {

        size_t inter_acc = std::accumulate(shots_per_time_.begin(), shots_per_time_.begin() + it, 0) * in_dimensions[0];
        in_dimensions.pop_back(); // Remove INT dimension
        in_dimensions.push_back(shots_per_time_[it]);

        for (size_t ich = 0; ich < CHA; ich++) {

            auto slice_view = cuNDArray<complext<REAL>>(in_dimensions, in->data() + stride_ch * ich + inter_acc);
            auto out_view_ch = cuNDArray<complext<REAL>>(out_dimensions, tmp.data() + stride_out * ich);

            plan_[it]->compute(slice_view, out_view_ch, &dcw_[it], NFFT_comp_mode::BACKWARDS_NC2C);
        }

        auto slice_view_output = cuNDArray<complext<REAL>>(out_dimensions, moving_images.data() + stride_out * it);

        this->mult_csm_conj_sum(&tmp, &slice_view_output);
        applyDeformation(&slice_view_output, forward_deformation_[it]);
        output += slice_view_output;
    }
    output /= complext<REAL>((REAL)shots_per_time_.size(),(REAL)shots_per_time_.size());
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::setup(_uint64d matrix_size, _uint64d matrix_size_os, REAL W) {
    for (auto ii = 0; ii < shots_per_time_.size(); ii++)
        plan_.push_back(NFFT<cuNDArray, REAL, D>::make_plan(matrix_size, matrix_size_os, W, convolutionType));
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::preprocess(std::vector<cuNDArray<_reald>>& trajectory) {
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
void cuNonCartesianMOCOOperator<REAL, D>::set_dcw(std::vector<cuNDArray<REAL>> dcw) {
    dcw_ = dcw;
}
template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::applyDeformation(cuNDArray<complext<REAL>> *moving_image, cuNDArray<REAL>  transformation) {

    // Setup solver
    boost::shared_ptr<cuLinearResampleOperator<REAL, D>> R(new cuLinearResampleOperator<REAL, D>());
    cuCKOpticalFlowSolver<REAL, D> CK;
    CK.set_interpolator(R);

    boost::shared_ptr<cuNDArray<REAL>> result = boost::make_shared<cuNDArray<REAL>>(transformation);

    auto mir = *real(moving_image);
    auto mii = *imag(moving_image);

    auto deformed_movingr = *CK.deform(&mir, result);
    auto deformed_movingi = *CK.deform(&mii, result);

    // Gadgetron::cuNDArray<complext<REAL>> deformed_image =
    //     *cureal_imag_to_complex<complext<REAL>>(&deformed_movingr, &deformed_movingi);
    moving_image = cureal_imag_to_complex<complext<REAL>>(&deformed_movingr, &deformed_movingi).get();
    //return deformed_image;
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::set_shots_per_time(std::vector<size_t> shots_per_time) {
    shots_per_time_ = shots_per_time;
}
template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::set_forward_deformation(std::vector<cuNDArray<REAL>> forward_deformation) {
    forward_deformation_ = forward_deformation;
}
template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::set_backward_deformation(std::vector<cuNDArray<REAL>> backward_deformation) {
    backward_deformation_ = backward_deformation;
}

template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 1>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 2>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 3>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 4>;

template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 1>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 2>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 3>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 4>;
