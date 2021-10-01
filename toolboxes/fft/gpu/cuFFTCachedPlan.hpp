#pragma once
#include "cuFFTCachedPlan.h"
#include "cuFFTPlan.h"

namespace {
namespace cufft_detail{
    auto compact_dims(int rank, const std::vector<size_t>& dimensions){
        assert(rank <= dimensions.size());
        auto dims = std::vector<size_t>(dimensions.begin(),dimensions.begin()+rank);
        auto batches = std::reduce(dimensions.begin()+rank,dimensions.end(),1,std::multiplies());
        dims.push_back(batches);
        return dims;
    }

    auto fetch_plan = [](auto& cache, const auto& compacted_dims ){
        if (!cache.count(compacted_dims)) cache.emplace(std::piecewise_construct,compacted_dims,{compacted_dims.size()-1,compacted_dims});
        return cache[compacted_dims];
    };
}
}

template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::fft1(cuNDArray<ComplexType>& in_out) {
    this->fft<1>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::fft2(cuNDArray<ComplexType>& in_out) {
    this->fft<2>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::fft3(cuNDArray<ComplexType>& in_out) {
    this->fft<3>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::ifft1(cuNDArray<ComplexType>& in_out) {
    this->ifft<1>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::ifft2(cuNDArray<ComplexType>& in_out) {
    this->ifft<2>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::ifft3(cuNDArray<ComplexType>& in_out) {
    this->ifft<3>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::fft1c(cuNDArray<ComplexType>& in_out) {
    this->fftc<1>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::fft2c(cuNDArray<ComplexType>& in_out) {
    this->fftc<2>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::fft3c(cuNDArray<ComplexType>& in_out) {
    this->fftc<3>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::ifft1c(cuNDArray<ComplexType>& in_out) {
    this->ifftc<1>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::ifft2c(cuNDArray<ComplexType>& in_out) {
    this->ifftc<2>(in_out);
}
template <class ComplexType> void Gadgetron::cuFFTCachedPlan<ComplexType>::ifft3c(cuNDArray<ComplexType>& in_out) {
    this->ifftc<3>(in_out);
}

template <class ComplexType>
template <int Rank>
void Gadgetron::cuFFTCachedPlan<ComplexType>::fft(Gadgetron::cuNDArray<ComplexType>& in_out) {
    cufft_detail::fetch_plan(plans,cufft_detail::compact_dims(Rank, in_out.dimensions())).fft(in_out);
}
template <class ComplexType>
template <int Rank>
void Gadgetron::cuFFTCachedPlan<ComplexType>::ifft(Gadgetron::cuNDArray<ComplexType>& in_out) {

    cufft_detail::fetch_plan(plans,cufft_detail::compact_dims(Rank, in_out.dimensions())).ifft(in_out);
}
template <class ComplexType>
template <int Rank>
void Gadgetron::cuFFTCachedPlan<ComplexType>::fftc(Gadgetron::cuNDArray<ComplexType>& in_out) {
    cufft_detail::fetch_plan(plans,cufft_detail::compact_dims(Rank, in_out.dimensions())).fftc(in_out);
}
template <class ComplexType>
template <int Rank>
void Gadgetron::cuFFTCachedPlan<ComplexType>::ifftc(Gadgetron::cuNDArray<ComplexType>& in_out) {
    cufft_detail::fetch_plan(plans,cufft_detail::compact_dims(Rank, in_out.dimensions())).ifftc(in_out);
}
