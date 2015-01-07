/** \file       gtPlusSPIRITNoNullSpace2DOperator.h
    \brief      Implement SPIRIT 2D operator without Null space
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRITNoNullSpaceOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRITNoNullSpace2DOperator : public gtPlusSPIRITNoNullSpaceOperator<T>
{
public:

    typedef gtPlusSPIRITNoNullSpaceOperator<T> BaseClass;

    gtPlusSPIRITNoNullSpace2DOperator() : BaseClass() {}
    virtual ~gtPlusSPIRITNoNullSpace2DOperator() {}

    virtual void printInfo(std::ostream& os);

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x);

protected:

};

template <typename T> 
void gtPlusSPIRITNoNullSpace2DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT 2D operator without null space constraint ------------------" << endl;
    os << "Implementation of SPIRIT 2D operator for ISMRMRD package" << endl;
    os << "-------------------------------------------------------------------------------------------------" << endl;
}

template <typename T> 
inline bool gtPlusSPIRITNoNullSpace2DOperator<T>::convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if ( !this->complexIm_Managed_.dimensions_equal(&x) )
    {
        this->complexIm_Managed_.create(x.get_dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, this->complexIm_Managed_);

    return true;
}

template <typename T> 
inline bool gtPlusSPIRITNoNullSpace2DOperator<T>::convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if ( !this->kspace_Managed_.dimensions_equal(&im) )
    {
        this->kspace_Managed_.create(im.get_dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, this->kspace_Managed_);

    return true;
}

}}
