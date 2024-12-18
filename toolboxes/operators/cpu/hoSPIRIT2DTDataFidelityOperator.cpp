
#include "hoSPIRIT2DTDataFidelityOperator.h"

namespace Gadgetron
{

template <typename T>
hoSPIRIT2DTDataFidelityOperator<T>::hoSPIRIT2DTDataFidelityOperator(const std::vector<size_t>& dims) : BaseClass(dims)
{
    this->use_non_centered_fft_ = false;
    this->no_null_space_ = true;
}

template <typename T>
hoSPIRIT2DTDataFidelityOperator<T>::~hoSPIRIT2DTDataFidelityOperator()
{
}

template <typename T>
void hoSPIRIT2DTDataFidelityOperator<T>::mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    if (!accumulate)
    {
        this->forward(*x, *y);
    }
    else
    {
        this->kspace_ = *y;
        this->forward(*x, *y);
        Gadgetron::add(this->kspace_, *y, *y);
    }
}

template <typename T>
void hoSPIRIT2DTDataFidelityOperator<T>::mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    if (!accumulate)
    {
        this->adjoint(*x, *y);
    }
    else
    {
        this->complexIm_ = *y;
        this->adjoint(*x, *y);
        Gadgetron::add(this->complexIm_, *y, *y);
    }
}

template <typename T>
void hoSPIRIT2DTDataFidelityOperator<T>::forward(const ARRAY_TYPE& x, ARRAY_TYPE& y)
{
    try
    {
        // y = [(G-I) + D]*x

        // x to image domain
        this->convert_to_image(x, this->complexIm_);

        // G-I
        this->apply_forward_kernel(this->complexIm_);

        this->convert_to_kspace(this->res_after_apply_kernel_sum_over_, y);

        // apply D, acquired points
        Gadgetron::multiply(this->acquired_points_indicator_, x, this->kspace_);

        Gadgetron::add(this->kspace_, y, y);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRIT2DTDataFidelityOperator<T>::forward(...) ... ");
    }
}

template <typename T>
void hoSPIRIT2DTDataFidelityOperator<T>::adjoint(const ARRAY_TYPE& x, ARRAY_TYPE& y)
{
    try
    {
        // y = [(G-I)' + D']*x

        // x to image domain
        this->convert_to_image(x, this->complexIm_dst_);

        // (G-I)'
        this->apply_adjoint_kernel(this->complexIm_dst_);

        this->convert_to_kspace(this->res_after_apply_kernel_sum_over_dst_, y);

        // apply D'
        Gadgetron::multiply(this->acquired_points_indicator_, x, this->kspace_);

        Gadgetron::add(this->kspace_, y, y);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRIT2DTDataFidelityOperator<T>::adjoint(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class hoSPIRIT2DTDataFidelityOperator< std::complex<float> >;
template class hoSPIRIT2DTDataFidelityOperator< std::complex<double> >;

}
