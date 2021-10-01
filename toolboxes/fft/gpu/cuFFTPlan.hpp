#pragma once
#include "cuNDFFT.h"
#include <numeric>

namespace Gadgetron::FFT_internal {
namespace {
template <class T> constexpr cufftType_t transform_type() {
    if constexpr (std::is_same_v<float, Gadgetron::realType_t<T>>)
        return CUFFT_C2C;
    if constexpr (std::is_same_v<double, Gadgetron::realType_t<T>>)
        return CUFFT_Z2Z;
}

std::string CUFFT_error_string(cufftResult error) {
    switch (error) {
    case CUFFT_SUCCESS:
        return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN:
        return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED:
        return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE:
        return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE:
        return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR:
        return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED:
        return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED:
        return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE:
        return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA:
        return "CUFFT_UNALIGNED_DATA";
    case CUFFT_INCOMPLETE_PARAMETER_LIST:
        return "CUFFT_INCOMPLETE_PARAMETER_LIST";
    case CUFFT_INVALID_DEVICE:
        return "CUFFT INVALID_DEVICE";
    case CUFFT_PARSE_ERROR:
        return "CUFFT_PARSE_ERROR";
    case CUFFT_NO_WORKSPACE:
        return "CUFFT_NO_WORKSPACE";
    case CUFFT_NOT_IMPLEMENTED:
        return "CUFFT_NOT_IMPLEMENTED";
    case CUFFT_LICENSE_ERROR:
        return "CUFFT_LICENSE_ERROR";
    case CUFFT_NOT_SUPPORTED:
        return "CUFFT_NOT_SUPPORTED";
    }

    return "<unknown>";
}

template <class ComplexType>
cufftResult_t executePlan(cufftHandle handle, const ComplexType* idata, ComplexType* odata, int direction) {
    if constexpr (std::is_same_v<float, Gadgetron::realType_t<ComplexType>>) {
        return cufftExecC2C(handle, (cufftComplex*)idata, (cufftComplex*)odata, direction);
    }
    if constexpr (std::is_same_v<double, Gadgetron::realType_t<ComplexType>>) {
        return cufftExecZ2Z(handle, (cufftDoubleComplex*)idata, (cufftDoubleComplex*)odata, direction);
    }
    return cufftResult::CUFFT_NOT_SUPPORTED;
}

template <class ComplexType> void timeswitch(cuNDArray<ComplexType>& in_out, int rank) {

    switch (rank) {
    case 1:
        return Gadgetron::timeswitch1D(&in_out);
    case 2:
        return Gadgetron::timeswitch2D(&in_out);
    default:
        Gadgetron::timeswitch3D(&in_out);
    }
    for (int i = 3; i < rank; i++)
        Gadgetron::timeswitch(&in_out, i);
}
} // namespace
}

template <class ComplexType, class ENABLER>
Gadgetron::cuFFTPlan<ComplexType, ENABLER>::cuFFTPlan(int rank, const std::vector<size_t>& dimensions) : rank(rank) {

    if (rank > dimensions.size())
        throw std::invalid_argument("Rank must be equal or smaller than the number of dimensions given ");

    auto int_dimensions = std::vector<int>(dimensions.begin(), dimensions.end());

    int dist = std::accumulate(dimensions.begin(), dimensions.end(), 1, std::multiplies());

    int batch_size = dimensions.size() > rank
                         ? std::accumulate(dimensions.begin() + rank, dimensions.end(), 1, std::multiplies())
                         : 1;

    auto result = cufftPlanMany(plan, rank, int_dimensions.data(), nullptr, 1, dist, int_dimensions.data(), 1, dist,
                                FFT_internal::transform_type<ComplexType>> (), batch_size);

    if (result != cufftResult::CUFFT_SUCCESS) {
        throw std::runtime_error("FFT plan failed with error" + FFT_internal::CUFFT_error_string(result));
    }
}
template <class ComplexType, class ENABLER> Gadgetron::cuFFTPlan<ComplexType, ENABLER>::~cuFFTPlan() {
    cufftDestroy(plan);
}

template <class ComplexType, class ENABLER>
void Gadgetron::cuFFTPlan<ComplexType, ENABLER>::fft(cuNDArray<ComplexType>& in_out) {
    if (!in_out.dimensions_equal(dimensions))
        throw std::runtime_error("Dimensions do not match FFT plan");
    auto result = FFT_internal::executePlan(plan, in_out.data(), in_out.data(), CUFFT_FORWARD);
    if (result != cufftResult::CUFFT_SUCCESS) {
        throw std::runtime_error("FFT failed with error" + FFT_internal::CUFFT_error_string(result));
    }
}

template <class ComplexType, class ENABLER>
void Gadgetron::cuFFTPlan<ComplexType, ENABLER>::ifft(cuNDArray<ComplexType>& in_out) {
    if (!in_out.dimensions_equal(dimensions))
        throw std::runtime_error("Dimensions do not match FFT plan");
    auto result = FFT_internal::executePlan(plan, in_out.data(), in_out.data(), CUFFT_INVERSE);
    if (result != cufftResult::CUFFT_SUCCESS) {
        throw std::runtime_error("IFFT failed with error" + FFT_internal::CUFFT_error_string(result));
    }
}


template <class ComplexType, class ENABLER>
void Gadgetron::cuFFTPlan<ComplexType, ENABLER>::fftc(cuNDArray<ComplexType>& in_out) {
    FFT_internal::timeswitch(in_out, rank);
    fft(in_out);
    FFT_internal::timeswitch(in_out, rank);
}
template <class ComplexType, class ENABLER>
void Gadgetron::cuFFTPlan<ComplexType, ENABLER>::ifftc(cuNDArray<ComplexType>& in_out) {
    FFT_internal::timeswitch(in_out, rank);
    ifft(in_out);
    FFT_internal::timeswitch(in_out, rank);
}
