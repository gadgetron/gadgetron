/** \file       gtPlusWaveletNoNullSpace3DOperator.h
    \brief      Implement 3D wavelet operator for without Null space cases
    \author     Hui Xue
*/

#pragma once

#include "gtPlusWavelet3DOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusWaveletNoNullSpace3DOperator : public gtPlusWavelet3DOperator<T>
{
public:

    typedef gtPlusWavelet3DOperator<T> BaseClass;

    gtPlusWaveletNoNullSpace3DOperator();
    virtual ~gtPlusWaveletNoNullSpace3DOperator();

    virtual void printInfo(std::ostream& os);

    // if the sensitivity S is set, compute gradient of ||wav*S'*F'*x||1
    // if not, compute gradient of ||wav*F'*x||1
    // x represents the unacquired kspace points [RO E1 CHA E2]
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // if the sensitivity S is set, compute cost value of L2 norm ||wav*S'*F'*x||1
    // if not, compute cost value of L2 norm ||wav*F'*x||1
    virtual bool obj(const hoNDArray<T>& x, T& obj);

    using BaseClass::scale_factor_first_dimension_;
    using BaseClass::scale_factor_second_dimension_;
    using BaseClass::scale_factor_third_dimension_;
    using BaseClass::change_coeffcients_third_dimension_boundary_;
    using BaseClass::numOfWavLevels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;
    using BaseClass::gtPlus_mem_manager_;

protected:

    // helper memory
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::res_after_apply_kernel_;
    using BaseClass::res_after_apply_kernel_sum_over_;

    using BaseClass::wav_coeff_norm_;
    using BaseClass::kspace_wav_;
    using BaseClass::complexIm_wav_;

    using BaseClass::kspace_Managed_;
    using BaseClass::complexIm_Managed_;
    using BaseClass::res_after_apply_kernel_Managed_;
    using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
gtPlusWaveletNoNullSpace3DOperator<T>::gtPlusWaveletNoNullSpace3DOperator() : BaseClass()
{

}

template <typename T> 
gtPlusWaveletNoNullSpace3DOperator<T>::~gtPlusWaveletNoNullSpace3DOperator()
{
}

template <typename T> 
inline bool gtPlusWaveletNoNullSpace3DOperator<T>::
grad(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(this->gradTask(x, g));
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletNoNullSpace3DOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
        return false;
    }

    return true;
}

template <typename T> 
inline bool gtPlusWaveletNoNullSpace3DOperator<T>::
obj(const hoNDArray<T>& x, T& obj)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(this->objTask(x, obj));
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletNoNullSpace3DOperator<T>::obj(const hoNDArray<T>& x, T& obj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
void gtPlusWaveletNoNullSpace3DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD wavelet 3D operator -----------------------" << endl;
    os << "Wavelet operator for gtPlus ISMRMRD package" << endl;
    os << "-------------------------------------------------------------------------" << endl;
}

}}
