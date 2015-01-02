/** \file       gtPlusWaveletNoNullSpace2DOperator.h
    \brief      Implement 2D wavelet operator for without Null space cases
    \author     Hui Xue
*/

#pragma once

#include "gtPlusWavelet2DOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusWaveletNoNullSpace2DOperator : public gtPlusWavelet2DOperator<T>
{
public:

    typedef gtPlusWavelet2DOperator<T> BaseClass;

    gtPlusWaveletNoNullSpace2DOperator();
    virtual ~gtPlusWaveletNoNullSpace2DOperator();

    virtual void printInfo(std::ostream& os);

    // if the sensitivity S is set, compute gradient of ||wav*S'*F'*x||1
    // if not, compute gradient of ||wav*F'*x||1
    // x represents the unacquired kspace points [RO E1 CHA]
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // if the sensitivity S is set, compute cost value of L2 norm ||wav*S'*F'*x||1
    // if not, compute cost value of L2 norm ||wav*F'*x||1
    virtual bool obj(const hoNDArray<T>& x, T& obj);

    using BaseClass::scale_factor_first_dimension_;
    using BaseClass::scale_factor_second_dimension_;
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
gtPlusWaveletNoNullSpace2DOperator<T>::gtPlusWaveletNoNullSpace2DOperator() : BaseClass()
{

}

template <typename T> 
gtPlusWaveletNoNullSpace2DOperator<T>::~gtPlusWaveletNoNullSpace2DOperator()
{
}

template <typename T> 
bool gtPlusWaveletNoNullSpace2DOperator<T>::
grad(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(this->gradTask(x, g));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletNoNullSpace2DOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWaveletNoNullSpace2DOperator<T>::
obj(const hoNDArray<T>& x, T& obj)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(this->objTask(x, obj));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletNoNullSpace2DOperator<T>::obj(const hoNDArray<T>& x, T& obj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
void gtPlusWaveletNoNullSpace2DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD wavelet 2D operator --------------------" << endl;
    os << "Wavelet 2D operator for gtPlus ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

}}
