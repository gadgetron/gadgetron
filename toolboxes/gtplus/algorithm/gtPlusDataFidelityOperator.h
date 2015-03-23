/** \file       gtPlusDataFidelityOperator.h
    \brief      Implement data fidelity operator
    \author     Hui Xue
*/

#pragma once

#include "gtPlusOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusDataFidelityOperator : public gtPlusOperator<T>
{
public:

    typedef gtPlusOperator<T> BaseClass;

    gtPlusDataFidelityOperator();
    virtual ~gtPlusDataFidelityOperator();

    virtual void printInfo(std::ostream& os);

    // forward operator
    // D*x
    virtual bool forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // adjoint operator
    // D'x
    virtual bool adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // gradient of ||Dx-y||2
    // 2*D'*(Dx-y)
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // L2 norm of ||Dx-y||2
    virtual bool obj(const hoNDArray<T>& x, T& obj);

    virtual bool unitary() const { return true; }

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;

protected:

    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::coil_senMap_;

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
gtPlusDataFidelityOperator<T>::gtPlusDataFidelityOperator() : BaseClass()
{

}

template <typename T> 
gtPlusDataFidelityOperator<T>::~gtPlusDataFidelityOperator()
{
}

template <typename T> 
bool gtPlusDataFidelityOperator<T>::
forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        Gadgetron::multiply(acquired_points_indicator_, x, y);
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusDataFidelityOperator<T>::forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusDataFidelityOperator<T>::
adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        Gadgetron::multiply(acquired_points_indicator_, x, y);
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusDataFidelityOperator<T>::adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusDataFidelityOperator<T>::
grad(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        // 2D'*(Dx-y)
        Gadgetron::multiply(acquired_points_indicator_, x, g);
        Gadgetron::subtract(g, *acquired_points_, g);
        Gadgetron::scal(T(2.0), g);
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusDataFidelityOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusDataFidelityOperator<T>::
obj(const hoNDArray<T>& x, T& obj)
{
    try
    {
        Gadgetron::multiply(acquired_points_indicator_, x, kspace_);
        Gadgetron::subtract(kspace_, *acquired_points_, kspace_);
        Gadgetron::dotc(kspace_, kspace_, obj);
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusDataFidelityOperator<T>::obj(const hoNDArray<T>& x, T& obj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
void gtPlusDataFidelityOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD data fidelity operator -----------------------" << endl;
    os << "Data fidelity operator for gtPlus ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------------" << endl;
}

}}
