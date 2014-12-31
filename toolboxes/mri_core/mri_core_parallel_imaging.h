
/** \file   mri_core_grappa.h
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{

template <typename T> 
class EXPORTMRICORE parallelImaging
{
public:

    typedef typename realType<T>::Type value_type;

    parallelImaging(){}
    virtual ~parallelImaging() {}

    virtual void printInfo(std::ostream& os);

    /// ---------------------------------------------------------------------
    /// 2D reconstruction
    /// ---------------------------------------------------------------------
    // compute unmixing coefficient from image domain kernel and coil sensitivity
    // kerIm: [RO E1 srcCHA dstCHA]
    // coilMap: [RO E1 dstCHA]
    // unmixCoeff: [RO E1 srcCHA]
    // gFactor: [RO E1]
    void unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, double acceFactorE1, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor);

    /// apply unmixing coefficient
    /// kspace: [RO E1 srcCHA ...]
    /// unmixCoeff : [RO E1 srcCHA]
    /// complexIm : [RO E1 ...]
    void applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);
    /// aliasedIm : [RO E1 srcCHA ...]
    void applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    /// ---------------------------------------------------------------------
    /// 3D reconstruction
    /// ---------------------------------------------------------------------
    // compute unmixing coefficient from image domain kernel and coil sensitivity
    // kerIm: [RO E1 E2 srcCHA dstCHA]
    // coilMap: [RO E1 E2 dstCHA]
    // unmixCoeff: [RO E1 E2 srcCHA]
    // gFactor: [RO E1 E2]
    void unmixCoeff3D(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, double acceFactorE1, double acceFactorE2, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor);

    /// apply unmixing coefficient
    /// kspace: [RO E1 E2 srcCHA ...]
    /// unmixCoeff : [RO E1 E2 srcCHA]
    /// complexIm : [RO E1 E2 ...]
    void applyUnmixCoeff3D(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);
    /// aliasedIm : [RO E1 E2 srcCHA ...]
    void applyUnmixCoeffImage3D(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

protected:

};

}
