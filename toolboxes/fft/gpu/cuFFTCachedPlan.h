#pragma once
#include "cuFFTPlan.h"
#include <map>

namespace Gadgetron {
template <class ComplexType> class cuFFTCachedPlan {
  public:
    cuFFTCachedPlan() = default;

    void fft1(cuNDArray<ComplexType>& in_out);
    void fft2(cuNDArray<ComplexType>& in_out);
    void fft3(cuNDArray<ComplexType>& in_out);

    void ifft1(cuNDArray<ComplexType>& in_out);
    void ifft2(cuNDArray<ComplexType>& in_out);
    void ifft3(cuNDArray<ComplexType>& in_out);

    void fft1c(cuNDArray<ComplexType>& in_out);
    void fft2c(cuNDArray<ComplexType>& in_out);
    void fft3c(cuNDArray<ComplexType>& in_out);

    void ifft1c(cuNDArray<ComplexType>& in_out);
    void ifft2c(cuNDArray<ComplexType>& in_out);
    void ifft3c(cuNDArray<ComplexType>& in_out);

    template <int Rank> void fft(cuNDArray<ComplexType>& in_out);
    template <int Rank> void ifft(cuNDArray<ComplexType>& in_out);

    template <int Rank> void fftc(cuNDArray<ComplexType>& in_out);

    template <int Rank> void ifftc(cuNDArray<ComplexType>& in_out);

  private:
    std::map<std::vector<size_t>, cuFFTPlan<ComplexType, void >> plans;
};

} // namespace Gadgetron
#include "cuFFTCachedPlan.hpp"