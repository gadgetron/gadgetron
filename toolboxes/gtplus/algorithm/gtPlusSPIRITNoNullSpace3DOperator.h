/** \file       gtPlusSPIRITNoNullSpace3DOperator.h
    \brief      Implement SPIRIT 3D operator without Null space
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRITNoNullSpaceOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRITNoNullSpace3DOperator : public gtPlusSPIRITNoNullSpaceOperator<T>
{
public:

    typedef gtPlusSPIRITNoNullSpaceOperator<T> BaseClass;

    gtPlusSPIRITNoNullSpace3DOperator() : BaseClass() {}
    virtual ~gtPlusSPIRITNoNullSpace3DOperator() {}

    virtual void printInfo(std::ostream& os);

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x);

protected:

};

template <typename T> 
void gtPlusSPIRITNoNullSpace3DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT 3D operator without null space constraint ------------------" << endl;
    os << "Implementation of SPIRIT 3D operator for ISMRMRD package" << endl;
    os << "-------------------------------------------------------------------------------------------------" << endl;
}

template <typename T> 
inline bool gtPlusSPIRITNoNullSpace3DOperator<T>::convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if ( !this->complexIm_Managed_.dimensions_equal(&x) )
    {
        this->complexIm_Managed_.create(x.get_dimensions());
    }

    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(x, im, this->complexIm_Managed_));
}

template <typename T> 
inline bool gtPlusSPIRITNoNullSpace3DOperator<T>::convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if ( !this->kspace_Managed_.dimensions_equal(&im) )
    {
        this->kspace_Managed_.create(im.get_dimensions());
    }

    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(im, x, this->kspace_Managed_));
}

}}
