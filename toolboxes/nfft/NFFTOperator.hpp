#include "NFFTOperator.h"
#include "complext.h"

namespace Gadgetron {



    template<template<class> class ARRAY, class REAL, unsigned int D>
    NFFTOperator<ARRAY, REAL, D>::NFFTOperator() : linearOperator<ARRAY < complext < REAL> > >() {
}


template<template<class> class ARRAY, class REAL, unsigned int D>
void
NFFTOperator<ARRAY, REAL, D>::mult_M(ARRAY <complext<REAL>> *in, ARRAY <complext<REAL>> *out, bool accumulate) {
    if (!in || !out) {
        throw std::runtime_error("NFFTOperator::mult_M : 0x0 input/output not accepted");
    }

    ARRAY <complext<REAL>> *tmp_out;

    if (accumulate) {
        tmp_out = new ARRAY <complext<REAL>>(out->dimensions());
    } else {
        tmp_out = out;
    }

    plan_->compute(*in, *tmp_out, dcw_.get(), NFFT_comp_mode::FORWARDS_C2NC);

    if (accumulate) {
        *out += *tmp_out;
        delete tmp_out;
    }
}

template<template<class> class ARRAY, class REAL, unsigned int D>
void
NFFTOperator<ARRAY, REAL, D>::mult_MH(ARRAY <complext<REAL>> *in, ARRAY <complext<REAL>> *out, bool accumulate) {
    if (!in || !out) {
        throw std::runtime_error("NFFTOperator::mult_MH : 0x0 input/output not accepted");
    }

    ARRAY <complext<REAL>> *tmp_out;

    if (accumulate) {
        tmp_out = new ARRAY <complext<REAL>>(out->dimensions());
    } else {
        tmp_out = out;
    }

    plan_->compute(*in, *tmp_out, dcw_.get(), NFFT_comp_mode::BACKWARDS_NC2C);
    if (accumulate) {
        *out += *tmp_out;
        delete tmp_out;
    }
}

template<template<class> class ARRAY, class REAL, unsigned int D>
void
NFFTOperator<ARRAY, REAL, D>::mult_MH_M(ARRAY <complext<REAL>> *in, ARRAY <complext<REAL>> *out, bool accumulate) {
    if (!in || !out) {
        throw std::runtime_error("NFFTOperator::mult_MH_M : 0x0 input/output not accepted");
    }

    boost::shared_ptr<std::vector<size_t> > codomain_dims = this->get_codomain_dimensions();
    if (codomain_dims.get() == 0x0 || codomain_dims->size() == 0) {
        throw std::runtime_error("NFFTOperator::mult_MH_M : operator codomain dimensions not set");
    }

    ARRAY <complext<REAL>> *tmp_out;

    if (accumulate) {
        tmp_out = new ARRAY <complext<REAL>>(out->dimensions());
    } else {
        tmp_out = out;
    }

    plan_->mult_MH_M(*in, *tmp_out, dcw_.get());

    if (accumulate) {
        *out += *tmp_out;
        delete tmp_out;
    }
}

template<template<class> class ARRAY, class REAL, unsigned int D>
void
NFFTOperator<ARRAY, REAL, D>::setup(typename uint64d<D>::Type matrix_size, typename uint64d<D>::Type matrix_size_os,
                               REAL W) {
    plan_ = NFFT<ARRAY,REAL,D>::make_plan(matrix_size, matrix_size_os, W);

}

template<template<class> class ARRAY, class REAL, unsigned int D>
void
NFFTOperator<ARRAY, REAL, D>::preprocess(const ARRAY<typename reald<REAL, D>::Type>& trajectory) {
    plan_->preprocess(trajectory, NFFT_prep_mode::ALL);
}

}
