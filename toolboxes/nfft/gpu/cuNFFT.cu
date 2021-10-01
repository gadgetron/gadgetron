/*
  CUDA implementation of the NFFT.

  -----------

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12): 1974-1985. 
*/

#include "cuNFFT.h"

#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_utils.h"
#include "cuNDFFT.h"

#include "NFFT.hpp"


using namespace Gadgetron;


template<class REAL, unsigned int D, ConvolutionType CONV>
Gadgetron::cuNFFT_impl<REAL, D, CONV>::cuNFFT_impl(
    const vector_td<size_t, D>& matrix_size,
    const vector_td<size_t, D>& matrix_size_os,
    REAL W,
    int device)
  : NFFT_plan<cuNDArray, REAL, D>(matrix_size, matrix_size_os, W)
{
    // Initialize gridding convolution. This was done in base class already, but
    // we need to do it again in order to provide the appropriate convolution
    // type.
    KaiserKernel<REAL, D> kernel(vector_td<unsigned int, D>(matrix_size),
                                 vector_td<unsigned int, D>(matrix_size_os),
                                 W);
    this->conv_ = GriddingConvolution<cuNDArray, complext<REAL>, D, KaiserKernel>::make(
        matrix_size, matrix_size_os, kernel, CONV);

    // Minimal initialization.
    this->initialize(device);
}


template<class REAL, unsigned int D, ConvolutionType CONV>
Gadgetron::cuNFFT_impl<REAL, D, CONV>::~cuNFFT_impl()
{

}


template<class REAL, unsigned int D, ConvolutionType CONV>
void
Gadgetron::cuNFFT_impl<REAL, D, CONV>::fft(cuNDArray <complext<REAL>>& data, NFFT_fft_mode mode, bool do_scale) {
    if (mode == NFFT_fft_mode::FORWARDS) {
        fft_plan.fft<D>(data);
    } else {
        fft_plan.ifft<D>(data);
    }

}


template<class REAL, unsigned int D, ConvolutionType CONV>
void
Gadgetron::cuNFFT_impl<REAL, D, CONV>::deapodize(cuNDArray <complext<REAL>>& image, bool fourier_domain) {

    typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>(*image.get_dimensions());
    bool oversampled_image = (image_dims == this->matrix_size_os_);

    if (!oversampled_image) {
        throw std::runtime_error("Error: cuNFFT_impl::deapodize: ERROR: oversampled image not provided as input.");
    }
    if (fourier_domain) {
        if (!deapodization_filterFFT)
            deapodization_filterFFT = compute_deapodization_filter(true);
        image *= *deapodization_filterFFT;
    } else {
        if (!deapodization_filter)
            deapodization_filter = compute_deapodization_filter(false);
        image *= *deapodization_filter;
    }
}


template<class REAL, unsigned int D, ConvolutionType CONV>
void Gadgetron::cuNFFT_impl<REAL, D, CONV>::initialize(int device)
{       
    // Device checks.
    if (cudaGetDevice(&this->device_) != cudaSuccess)
    {
        throw cuda_error("Error: cuNFFT_impl::barebones:: unable to get this->device_ no");
    }
    
    if (device < 0)
    {
        if (cudaGetDevice(&this->device_) != cudaSuccess)
        {
            throw cuda_error("Error: cuNFFT_impl::setup: unable to determine "
                             "device properties.");
        }
    }
    else
    {
        this->device_ = device;
    }

    int device_no_old;
    if (cudaGetDevice(&device_no_old) != cudaSuccess)
        throw cuda_error("Error: cuNFFT_impl::setup: unable to get device number");

    if (this->device_ != device_no_old &&
        cudaSetDevice(this->device_) != cudaSuccess)
        throw cuda_error("Error: cuNFFT_impl::setup: unable to set device");

    if (this->device_ != device_no_old &&
        cudaSetDevice(device_no_old) != cudaSuccess)
        throw cuda_error("Error: cuNFFT_impl::setup: unable to restore device");
}


//
// Grid fictitious trajectory with a single sample at the origin
//

template<class REAL, unsigned int D, template<class, unsigned int> class K>
__global__ void
compute_deapodization_filter_kernel(typename uintd<D>::Type matrix_size_os,
                                    typename reald<REAL, D>::Type matrix_size_os_real,
                                    complext <REAL> *__restrict__ image_os,
                                    const ConvolutionKernel<REAL, D, K>* kernel)
{
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int num_elements = prod(matrix_size_os);

    if (idx < num_elements) {

        // Compute weight from Kaiser-Bessel filter
        const typename uintd<D>::Type cell_pos = idx_to_co(idx, matrix_size_os);

        // Sample position ("origin")
        const vector_td<REAL, D> sample_pos = REAL(0.5) * matrix_size_os_real;

        // Calculate the distance between the cell and the sample
        vector_td<REAL, D>
        cell_pos_real = vector_td<REAL, D>(cell_pos);
        const typename reald<REAL, D>::Type delta = abs(sample_pos - cell_pos_real);

        // Compute convolution weight.
        REAL weight = kernel->get(delta);

        // Output weight
        image_os[idx] = complext<REAL>(weight, 0.0f);
    }
}

//
// Function to calculate the deapodization filter
//

template<class REAL, unsigned int D, ConvolutionType CONV>
boost::shared_ptr<cuNDArray < complext < REAL> > >
Gadgetron::cuNFFT_impl<REAL, D, CONV>::compute_deapodization_filter(bool FFTed) {
    std::vector<size_t> tmp_vec_os = to_std_vector(this->matrix_size_os_);

    auto  filter = boost::make_shared<cuNDArray<complext<REAL>>>(tmp_vec_os);
    vector_td<REAL, D>
    matrix_size_os_real = vector_td<REAL, D>(this->matrix_size_os_);

    // Find dimensions of grid/blocks.
    dim3 dimBlock(256);
    dim3 dimGrid((prod(this->matrix_size_os_) + dimBlock.x - 1) / dimBlock.x);

    // Invoke kernel
    compute_deapodization_filter_kernel<REAL, D><<<dimGrid, dimBlock>>>(
        vector_td<unsigned int, D>(
        this->matrix_size_os_), matrix_size_os_real,
        filter->get_data_ptr(),
        dynamic_cast<const cuGriddingConvolution<complext<REAL>, D, KaiserKernel>*>(
            this->conv_.get())->get_kernel_d());

    CHECK_FOR_CUDA_ERROR();

    // FFT
    if (FFTed) {
        fft(*filter, NFFT_fft_mode::FORWARDS, false);
    } else {
        fft(*filter, NFFT_fft_mode::BACKWARDS, false);
    }
    // Reciprocal
    reciprocal_inplace(filter.get());
    return filter;
}


namespace Gadgetron {
    template<class REAL, unsigned int D>
    boost::shared_ptr<cuNFFT_plan<REAL, D>> NFFT<cuNDArray, REAL, D>::make_plan(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os, REAL W, ConvolutionType conv) {
        switch (conv) {
            case ConvolutionType::STANDARD:
                return boost::make_shared<cuNFFT_impl<REAL, D, ConvolutionType::STANDARD> >(matrix_size,matrix_size_os,W);
            case ConvolutionType::ATOMIC:
                return boost::make_shared<cuNFFT_impl<REAL, D, ConvolutionType::ATOMIC> >(matrix_size,matrix_size_os,W);
            case ConvolutionType::SPARSE_MATRIX:
                return boost::make_shared<cuNFFT_impl<REAL, D, ConvolutionType::SPARSE_MATRIX>>(matrix_size,matrix_size_os,W);
        }
        throw std::runtime_error(
                "Invalid convolution type provided. If you're reading this, you may have broken your computer quite badly");
    }

    template<unsigned int D>

    boost::shared_ptr<cuNFFT_plan<double, D>> NFFT<cuNDArray, double, D>::make_plan(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os, double W, ConvolutionType conv) {
        if (conv == ConvolutionType::STANDARD) {
            return boost::make_shared<cuNFFT_impl<double, D, ConvolutionType::STANDARD>>(matrix_size,matrix_size_os,W);
        }
        throw std::runtime_error("Only standard convolution type supported for doubles");
    }
}

//
// Template instantion
//

template
class Gadgetron::cuNFFT_impl<float, 1, ConvolutionType::ATOMIC>;


template
class Gadgetron::cuNFFT_impl<float, 1, ConvolutionType::SPARSE_MATRIX>;

template
class Gadgetron::cuNFFT_impl<float, 1>;

template
class Gadgetron::cuNFFT_impl<double, 1>;

template
class Gadgetron::cuNFFT_impl<float, 2, ConvolutionType::ATOMIC>;


template
class Gadgetron::cuNFFT_impl<float, 2, ConvolutionType::SPARSE_MATRIX>;

template
class Gadgetron::cuNFFT_impl<float, 2>;

template
class Gadgetron::cuNFFT_impl<double, 2>;

template
class Gadgetron::cuNFFT_impl<float, 3, ConvolutionType::ATOMIC>;

template
class Gadgetron::cuNFFT_impl<float, 3, ConvolutionType::SPARSE_MATRIX>;

template
class Gadgetron::cuNFFT_impl<float, 3>;

template
class Gadgetron::cuNFFT_impl<double, 3>;

template
class Gadgetron::cuNFFT_impl<float, 4, ConvolutionType::ATOMIC>;


template
class Gadgetron::cuNFFT_impl<float, 4, ConvolutionType::SPARSE_MATRIX>;
template
class Gadgetron::cuNFFT_impl<float, 4>;

template
class Gadgetron::cuNFFT_impl<double, 4>;


template class Gadgetron::NFFT<cuNDArray,float,1>;
template class Gadgetron::NFFT<cuNDArray,float,2>;
template class Gadgetron::NFFT<cuNDArray,float,3>;
template class Gadgetron::NFFT<cuNDArray,float,4>;



template class Gadgetron::NFFT<cuNDArray,double,1>;
template class Gadgetron::NFFT<cuNDArray,double,2>;
template class Gadgetron::NFFT<cuNDArray,double,3>;
template class Gadgetron::NFFT<cuNDArray,double,4>;
