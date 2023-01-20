#pragma once
#include <cufft.h>
#include "complext.h"
namespace Gadgetron {

template <class ComplexType, class = std::enable_if_t<is_complex_type_v<ComplexType>>> class cuFFTPlan {

  public:
    /**
     *
     * @param rank Dimensionality of the FFT, i.e 1, 2 or 3
     * @param dimensions Domain size of the FFT. Must be at least of length rank. Will batch over further dimensions
     */
    cuFFTPlan(int rank, const std::vector<size_t>& dimensions);

    ~cuFFTPlan();

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
    
    /**
     * Creates a non-centered inplace FFT
     * @param in_out
     */
    void fft(cuNDArray<ComplexType>& in_out, bool scale = true);
;
    /**
    * Creates a non-centered inplace inverse FFT
    * @param in_out
     */
    void ifft(cuNDArray<ComplexType>& in_out, bool scale = true);
;

    /**
    * Creates a centered inplace FFT
    * @param in_out
     */
    void fftc(cuNDArray<ComplexType>& in_out, bool scale=true );

    /**
    * Created a centered inplace inverse FFT
    * @param in_out
     */
    void ifftc(cuNDArray<ComplexType>& in_out, bool scale=true);
;

  private:
    cufftHandle plan;
    const int rank;
    const std::vector<size_t> dimensions;
};
}

#include "cuFFTPlan.hpp"
