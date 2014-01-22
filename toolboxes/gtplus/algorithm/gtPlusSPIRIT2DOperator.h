/** \file       gtPlusSPIRIT2DOperator.h
    \brief      Base class for gtPlus 2D operators
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRITOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRIT2DOperator : public gtPlusSPIRITOperator<T>
{
public:

    typedef gtPlusSPIRITOperator<T> BaseClass;

    gtPlusSPIRIT2DOperator() : BaseClass() {}
    virtual ~gtPlusSPIRIT2DOperator() {}

    virtual void printInfo(std::ostream& os);

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x);

    using BaseClass::use_symmetric_spirit_;
    using BaseClass::use_non_centered_fft_;
    using BaseClass::calib_use_gpu_;

protected:

    // [RO E1 srcCHA dstCHA]
    using BaseClass::forward_kernel_;
    using BaseClass::adjoint_kernel_;
    using BaseClass::adjoint_forward_kernel_;
    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;

    // helper memory
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::res_after_apply_kernel_;
    using BaseClass::res_after_apply_kernel_sum_over_;

    using BaseClass::kspace_Managed_;
    using BaseClass::complexIm_Managed_;
    using BaseClass::res_after_apply_kernel_Managed_;
    using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
void gtPlusSPIRIT2DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT 2D operator ------------------" << endl;
    os << "Implementation of SPIRIT 2D operator for ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
inline bool gtPlusSPIRIT2DOperator<T>::convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if ( this->use_non_centered_fft_ )
    {
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2(x, im));
    }
    else
    {
        if ( !complexIm_Managed_.dimensions_equal(&x) )
        {
            complexIm_Managed_.create(x.get_dimensions());
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, complexIm_Managed_));
    }

    return true;
}

template <typename T> 
inline bool gtPlusSPIRIT2DOperator<T>::convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if ( this->use_non_centered_fft_ )
    {
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2(im, x));
    }
    else
    {
        if ( !kspace_Managed_.dimensions_equal(&im) )
        {
            kspace_Managed_.create(im.get_dimensions());
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, kspace_Managed_));
    }

    return true;
}

}}
