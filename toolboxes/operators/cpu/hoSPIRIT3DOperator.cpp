
#include "hoSPIRIT3DOperator.h"
#include "hoNDFFT.h"

namespace Gadgetron 
{

template <typename T> 
hoSPIRIT3DOperator<T>::hoSPIRIT3DOperator(std::vector<size_t> *dims) : BaseClass(dims)
{
}

template <typename T> 
hoSPIRIT3DOperator<T>::~hoSPIRIT3DOperator()
{
}

template <typename T>
void hoSPIRIT3DOperator<T>::convert_to_image(const ARRAY_TYPE& x, ARRAY_TYPE& im)
{
    try
    {
        if (this->use_non_centered_fft_)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3(x, im);
        }
        else
        {
            if (!fft_im_buffer_.dimensions_equal(x))
            {
                fft_im_buffer_.create(x.dimensions());
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(x, im, fft_im_buffer_);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in hoSPIRIT3DOperator<T>::convert_to_image(...) ... ");
    }
}

template <typename T>
void hoSPIRIT3DOperator<T>::convert_to_kspace(const ARRAY_TYPE& im, ARRAY_TYPE& x)
{
    try
    {
        if (this->use_non_centered_fft_)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3(im, x);
        }
        else
        {
            if (!fft_kspace_buffer_.dimensions_equal(im))
            {
                fft_kspace_buffer_.create(im.dimensions());
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(im, x, fft_kspace_buffer_);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in hoSPIRIT3DOperator<T>::convert_to_kspace(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCPUOPERATOR hoSPIRIT3DOperator< std::complex<float> >;
template class EXPORTCPUOPERATOR hoSPIRIT3DOperator< std::complex<double> >;

}
