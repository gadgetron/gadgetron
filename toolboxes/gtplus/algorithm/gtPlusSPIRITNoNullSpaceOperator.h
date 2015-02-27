/** \file       gtPlusSPIRITNoNullSpaceOperator.h
    \brief      Base class for SPIRIT operators without Null space
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRITOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRITNoNullSpaceOperator : public gtPlusSPIRITOperator<T>
{
public:

    typedef gtPlusSPIRITOperator<T> BaseClass;

    gtPlusSPIRITNoNullSpaceOperator() : BaseClass() {}
    virtual ~gtPlusSPIRITNoNullSpaceOperator() {}

    virtual void printInfo(std::ostream& os);

    // compute gradient of ||(G-I)x||2
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // compute cost value of L2 norm ||(G-I)x||2
    virtual bool obj(const hoNDArray<T>& x, T& obj);

    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::performTiming_;
    //using BaseClass::gt_exporter_;
    //using BaseClass::debugFolder_;
    //using BaseClass::gtPlus_util_;
    //using BaseClass::gtPlus_util_complex_;
    //using BaseClass::gtPlus_mem_manager_;
    //using BaseClass::use_symmetric_spirit_;

protected:

    // [... srcCHA dstCHA]
    //using BaseClass::forward_kernel_;
    //using BaseClass::adjoint_kernel_;
    //using BaseClass::adjoint_forward_kernel_;
    //using BaseClass::acquired_points_;
    //using BaseClass::acquired_points_indicator_;
    //using BaseClass::unacquired_points_indicator_;

    // helper memory
    //using BaseClass::kspace_;
    //using BaseClass::complexIm_;
    //using BaseClass::res_after_apply_kernel_;
    //using BaseClass::res_after_apply_kernel_sum_over_;

    //using BaseClass::kspace_Managed_;
    //using BaseClass::complexIm_Managed_;
    //using BaseClass::res_after_apply_kernel_Managed_;
    //using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
void gtPlusSPIRITNoNullSpaceOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT operator without null space constraint ------------------" << endl;
    os << "Implementation of SPIRIT operator for ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusSPIRITNoNullSpaceOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        // gradient of L2 norm is
        // 2*(G-I)'(G-I)x

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, this->complexIm_));

        // apply kernel and sum
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*this->adjoint_forward_kernel_, this->complexIm_, this->res_after_apply_kernel_));
        GADGET_CHECK_RETURN_FALSE(this->performSumOverSrcChannel(this->res_after_apply_kernel_, this->res_after_apply_kernel_sum_over_));

        // go back to kspace 
        GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(this->res_after_apply_kernel_sum_over_, g));

        // multiply by 2
        Gadgetron::scal(T(2.0), g);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITNoNullSpaceOperator<T>::grad(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITNoNullSpaceOperator<T>::obj(const hoNDArray<T>& x, T& obj)
{
    try
    {
        // L2 norm
        // ||(G-I)x||2

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, this->complexIm_));

        // apply kernel and sum
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*this->forward_kernel_, this->complexIm_, this->res_after_apply_kernel_));
        GADGET_CHECK_RETURN_FALSE(this->performSumOverSrcChannel(this->res_after_apply_kernel_, this->res_after_apply_kernel_sum_over_));

        // L2 norm
        Gadgetron::dotc(this->res_after_apply_kernel_sum_over_, this->res_after_apply_kernel_sum_over_, obj);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITNoNullSpaceOperator<T>::grad(...) ... ");
        return false;
    }

    return true;
}

}}
