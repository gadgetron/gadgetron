/** \file       gtPlusSPIRIT3DOperator.h
    \brief      Base class for gtPlus 3D operators
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRITOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRIT3DOperator : public gtPlusSPIRITOperator<T>
{
public:

    typedef gtPlusSPIRITOperator<T> BaseClass;

    gtPlusSPIRIT3DOperator() : BaseClass() {}
    virtual ~gtPlusSPIRIT3DOperator() {}

    virtual void printInfo(std::ostream& os);

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x);

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;
    using BaseClass::gtPlus_mem_manager_;
    using BaseClass::use_symmetric_spirit_;

protected:

    // [RO E1 srcCHA dstCHA]
    using BaseClass::forward_kernel_;
    using BaseClass::djoint_kernel_;
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
void gtPlusSPIRIT3DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT 3D operator ------------------" << endl;
    os << "Implementation of SPIRIT 3D operator for ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
inline bool gtPlusSPIRIT3DOperator<T>::convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if ( !complexIm_Managed_.dimensions_equal(&x) )
    {
        complexIm_Managed_.create(x.get_dimensions());
    }

    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(x, im, complexIm_Managed_));

    return true;
}

template <typename T> 
inline bool gtPlusSPIRIT3DOperator<T>::convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if ( !kspace_Managed_.dimensions_equal(&im) )
    {
        kspace_Managed_.create(im.get_dimensions());
    }

    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(im, x, kspace_Managed_));

    return true;
}

}}
