#pragma once
#include "cuFFTPlan.h"
#include <map>

namespace Gadgetron {
template <class ComplexType> class cuFFTCachedPlan {
  public:
    cuFFTCachedPlan() = default;

    void fft1(cuNDArray<ComplexType>& in_out, bool scale=true);
    void fft2(cuNDArray<ComplexType>& in_out, bool scale=true);
    void fft3(cuNDArray<ComplexType>& in_out, bool scale=true);

    void ifft1(cuNDArray<ComplexType>& in_out, bool scale=true);
    void ifft2(cuNDArray<ComplexType>& in_out, bool scale=true);
    void ifft3(cuNDArray<ComplexType>& in_out, bool scale=true);

    void fft1c(cuNDArray<ComplexType>& in_out, bool scale=true);
    void fft2c(cuNDArray<ComplexType>& in_out, bool scale=true);
    void fft3c(cuNDArray<ComplexType>& in_out, bool scale=true);

    void ifft1c(cuNDArray<ComplexType>& in_out, bool scale=true);
    void ifft2c(cuNDArray<ComplexType>& in_out, bool scale=true);
    void ifft3c(cuNDArray<ComplexType>& in_out, bool scale=true);

    void fft(cuNDArray<ComplexType>& in_out, unsigned int rank, bool scale=true);
    void ifft(cuNDArray<ComplexType>& in_out, unsigned int rank, bool scale=true);

    void fftc(cuNDArray<ComplexType>& in_out, unsigned int rank, bool scale=true);
    void ifftc(cuNDArray<ComplexType>& in_out, unsigned int rank, bool scale=true);

  private:
    std::map<std::vector<size_t>, cuFFTPlan<ComplexType, void >> plans;
};

} // namespace Gadgetron
#include "cuFFTCachedPlan.hpp"