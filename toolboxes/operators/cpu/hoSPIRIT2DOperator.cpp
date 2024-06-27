
#include "hoSPIRIT2DOperator.h"
#include "hoNDFFT.h"

namespace Gadgetron 
{

template <typename T> 
hoSPIRIT2DOperator<T>::hoSPIRIT2DOperator(std::vector<size_t> *dims) : BaseClass(dims)
{
}

template <typename T> 
hoSPIRIT2DOperator<T>::~hoSPIRIT2DOperator()
{
}

template <typename T>
void hoSPIRIT2DOperator<T>::convert_to_image(const ARRAY_TYPE& x, ARRAY_TYPE& im)
{
    try
    {
        if (this->use_non_centered_fft_)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2(x, im);
        }
        else
        {
            if (!fft_im_buffer_.dimensions_equal(x))
            {
                fft_im_buffer_.create(x.dimensions());
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, fft_im_buffer_);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in hoSPIRIT2DOperator<T>::convert_to_image(...) ... ");
    }
}

template <typename T>
void hoSPIRIT2DOperator<T>::convert_to_kspace(const ARRAY_TYPE& im, ARRAY_TYPE& x)
{
    try
    {
        if (this->use_non_centered_fft_)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2(im, x);
        }
        else
        {
            if (!fft_kspace_buffer_.dimensions_equal(im))
            {
                fft_kspace_buffer_.create(im.dimensions());
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, fft_kspace_buffer_);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in hoSPIRIT2DOperator<T>::convert_to_kspace(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCPUOPERATOR hoSPIRIT2DOperator< std::complex<float> >;
template class EXPORTCPUOPERATOR hoSPIRIT2DOperator< std::complex<double> >;

}
