/** \file       hoSPIRIT2DTDataFidelityOperator.h
    \brief      SPIRIT 2DT operator with data fidelity
    \author     Hui Xue
*/

#pragma once

#include "hoSPIRIT2DTOperator.h"

namespace Gadgetron { 

template <typename T> 
class EXPORTCPUOPERATOR hoSPIRIT2DTDataFidelityOperator : public hoSPIRIT2DTOperator<T>
{
public:

    typedef hoSPIRIT2DTOperator<T> BaseClass;
    typedef typename BaseClass::ARRAY_TYPE ARRAY_TYPE;
    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;

    hoSPIRIT2DTDataFidelityOperator(const std::vector<size_t>& dims);
    virtual ~hoSPIRIT2DTDataFidelityOperator();

    // x is the kspace, including acquired and unacquired points
    // y = [(G-I) + D]*x
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

    // y = [(G-I)' + D']*x
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

protected:

    void forward(const ARRAY_TYPE& x, ARRAY_TYPE& y);
    void adjoint(const ARRAY_TYPE& x, ARRAY_TYPE& y);
};

}
