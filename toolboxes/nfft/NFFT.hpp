#include "NFFT.h"

namespace Gadgetron {

    template<template<class> class ARRAY, class REAL, unsigned int D>
    NFFT_plan<ARRAY, REAL, D>::NFFT_plan(const vector_td <size_t, D> &matrix_size,
                                         const vector_td <size_t, D> &matrix_size_os, REAL W) : matrix_size(
            matrix_size), matrix_size_os(matrix_size_os), W(W) {
        if (W < REAL(1.0))
            throw std::runtime_error("Kernel width must be larger than 1");
        if (matrix_size > matrix_size_os)
            throw std::runtime_error("Oversampled matrix size must be as least as great as the matrix size");
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    NFFT_plan<ARRAY, REAL, D>::NFFT_plan(const vector_td <size_t, D> &matrix_size, REAL oversampling_factor, REAL W)
            : matrix_size(matrix_size), W(W) {

        matrix_size_os = vector_td<size_t, D>(vector_td<REAL, D>(matrix_size) * oversampling_factor);
        if (W < REAL(1.0))
            throw std::runtime_error("Kernel width must be larger than 1");
        if (oversampling_factor <= REAL(1.0))
            throw std::runtime_error("Oversampling factor must be 1 or larger");
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    void NFFT_plan<ARRAY, REAL, D>::preprocess(const ARRAY<vector_td<REAL,D>> &trajectory, NFFT_prep_mode mode) {
        number_of_samples = trajectory.get_size(0);
        number_of_frames = trajectory.get_number_of_elements()/number_of_samples;

    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    void NFFT_plan<ARRAY, REAL, D>::mult_MH_M(const ARRAY<complext<REAL>>& in, ARRAY<complext<REAL>>& out, const ARRAY<REAL>* dcw) {

        size_t number_of_batches = in.get_number_of_elements()/(prod(this->matrix_size)*this->number_of_frames);
        ARRAY<complext<REAL>> tmp(this->number_of_samples,this->number_of_frames,number_of_batches);
        compute(in, tmp, dcw, NFFT_comp_mode::FORWARDS_C2NC);
        compute(tmp, out,dcw, NFFT_comp_mode::BACKWARDS_NC2C);

    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
void NFFT_plan<ARRAY,REAL,D>::compute(const ARRAY<complext<REAL>>& in, ARRAY<complext<REAL>>&out,
        const ARRAY<REAL> *dcw, NFFT_comp_mode mode){

    typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>
            ((mode == NFFT_comp_mode::FORWARDS_C2NC || mode == NFFT_comp_mode::BACKWARDS_C2NC) ? *in.get_dimensions()
                                                                                               : *out.get_dimensions());
    bool oversampled_image = (image_dims == this->matrix_size_os);

    auto vec_dims = to_std_vector(this->matrix_size_os);
    {
        const auto image = ((mode == NFFT_comp_mode::FORWARDS_C2NC ||
                             mode == NFFT_comp_mode::BACKWARDS_C2NC) ? &in : &out);
        for (unsigned int d = D; d<image->get_number_of_dimensions(); d++)
            vec_dims.push_back(image->get_size(d));
    }

    switch (mode) {

        case NFFT_comp_mode::FORWARDS_C2NC: {
            if (!oversampled_image) {
                auto working_image = ARRAY<complext<REAL >> (vec_dims);
                pad<complext<REAL >,D> (in, working_image);
                compute_NFFT_C2NC(working_image, out);
            } else {
                auto copy = in;
                compute_NFFT_C2NC(copy, out);
            }

            if (dcw) {
                out *= *dcw;
            }

        }
            break;

        case NFFT_comp_mode::FORWARDS_NC2C: {

            // Density compensation
            auto working_samples = in;
            boost::shared_ptr<ARRAY<complext<REAL >> > samples_dcw;
            if (dcw) {
                working_samples *= *dcw;
            }


            if (!oversampled_image) {
                auto working_image = ARRAY<complext<REAL >> (&vec_dims);
                compute_NFFT_NC2C(working_samples, working_image);
                crop<complext<REAL > , D >
                                           ((this->matrix_size_os - this->matrix_size) >> 1,this->matrix_size, working_image, out);
            } else {
                compute_NFFT_NC2C(working_samples, out);
            }
        }
            break;

        case NFFT_comp_mode::BACKWARDS_NC2C: {

            // Density compensation
            const ARRAY<complext<REAL>> *working_samples = &in;
            boost::shared_ptr<ARRAY<complext<REAL >> > samples_dcw;
            if (dcw) {
                samples_dcw = boost::make_shared<ARRAY<complext<REAL >> > (in);
                *samples_dcw *= *dcw;
                working_samples = samples_dcw.get();
            }


            if (!oversampled_image) {
                auto working_image = ARRAY<complext<REAL>>(vec_dims);
                compute_NFFTH_NC2C(*working_samples, working_image);
                crop<complext<REAL >,D>((this->matrix_size_os - this->matrix_size) >> 1,this->matrix_size, working_image, out);
            } else {
                compute_NFFTH_NC2C(*working_samples, out);
            }

        }
            break;

        case NFFT_comp_mode::BACKWARDS_C2NC: {

            if (!oversampled_image) {
                auto working_image = ARRAY<complext<REAL >> (vec_dims);
                pad<complext<REAL >,D>(in, working_image);
                compute_NFFTH_C2NC(working_image, out);
            } else {
                auto copy = in;
                compute_NFFTH_C2NC(copy, out);
            }

            if (dcw)
                out *= *dcw;
        }

            break;
    };
}


template<template<class> class ARRAY, class REAL, unsigned int D>
void NFFT_plan<ARRAY, REAL, D>::compute_NFFT_C2NC(ARRAY<complext<REAL>>& image, ARRAY<complext<REAL>>& samples) {

    deapodize(image);
    fft(image, NFFT_fft_mode::FORWARDS);
    convolve(image, samples,  NFFT_conv_mode::C2NC);
}

template<template<class> class ARRAY, class REAL, unsigned int D>
void NFFT_plan<ARRAY, REAL, D>::compute_NFFT_NC2C(const ARRAY<complext<REAL>>& samples, ARRAY<complext<REAL>>& image) {

    convolve(samples, image, NFFT_conv_mode::NC2C);
    fft(image, NFFT_fft_mode::FORWARDS);
    deapodize(image, true);
}

template<template<class> class ARRAY, class REAL, unsigned int D>
void NFFT_plan<ARRAY, REAL, D>::compute_NFFTH_NC2C(const ARRAY<complext<REAL>>& samples, ARRAY<complext<REAL>>& image) {

    convolve(samples, image,  NFFT_conv_mode::NC2C);
    fft(image, NFFT_fft_mode::BACKWARDS);
    deapodize(image);
}

template<template<class> class ARRAY, class REAL, unsigned int D>
void NFFT_plan<ARRAY, REAL, D>::compute_NFFTH_C2NC(ARRAY<complext<REAL>>& image, ARRAY<complext<REAL>>& samples) {

    deapodize(image, true);
    convolve(samples, image, NFFT_conv_mode::NC2C);
    fft(image, NFFT_fft_mode::BACKWARDS);
}


}
