
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
void hoSPIRIT2DOperator<T>::convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    try
    {
        if (this->use_non_centered_fft_)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2(x, im);
        }
        else
        {
            if (!complexIm_.dimensions_equal(&x))
            {
                complexIm_.create(x.get_dimensions());
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, complexIm_);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in hoSPIRIT2DOperator<T>::convert_to_image(...) ... ");
    }
}

template <typename T>
void hoSPIRIT2DOperator<T>::convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    try
    {
        if (this->use_non_centered_fft_)
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2(im, x);
        }
        else
        {
            if (!kspace_.dimensions_equal(&im))
            {
                kspace_.create(im.get_dimensions());
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, kspace_);
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
