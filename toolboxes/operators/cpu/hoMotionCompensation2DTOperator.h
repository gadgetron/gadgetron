/** \file   hoMotionCompensation2DTOperator.h
    \brief  Implement motion compensation operator for 2DT cases
    \author Hui Xue
*/

#pragma once

#include "hoNDArray_linalg.h"
#include "hoImageRegContainer2DRegistration.h"

namespace Gadgetron {

template <typename T, typename CoordType>
class hoMotionCompensation2DTOperator
{
public:

    typedef typename realType<T>::Type value_type;

    typedef hoNDArray<T> ARRAY_TYPE;

    hoMotionCompensation2DTOperator();
    virtual ~hoMotionCompensation2DTOperator();

    /// x: [RO E1 CHA N]
    /// y: [RO E1 CHA N]
    /// apply forward deformation Mx
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);
    /// appy adjoint deformation M'x
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

    /// compute gradient of ||Mx||1
    virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate = false);

    /// compute cost value of L1 norm ||Mx||1
    virtual value_type magnitude(ARRAY_TYPE* x);

    T bg_value_;
    Gadgetron::GT_BOUNDARY_CONDITION bc_;

    hoNDArray<CoordType> dx_;
    hoNDArray<CoordType> dy_;
    hoNDArray<CoordType> adj_dx_;
    hoNDArray<CoordType> adj_dy_;

protected:

    // warp the 2D+T image arrays
    // im : complex image [RO E1 CHA N]
    // dx, dy : [RO E1 N] deformation fields
    // the image domain warpping is performed
    virtual void warp_image(const hoNDArray<T>& im, const hoNDArray<CoordType>& dx, const hoNDArray<CoordType>& dy, hoNDArray<T>& warpped, T bgValue=0, Gadgetron::GT_BOUNDARY_CONDITION bh=GT_BOUNDARY_CONDITION_FIXEDVALUE);

    // helper memory
    ARRAY_TYPE moco_im_;
    ARRAY_TYPE adj_moco_im_;

    ARRAY_TYPE grad_im_;

    hoNDArray<value_type> moco_im_Norm_;
    hoNDArray<value_type> moco_im_Norm_approx_;
};

}
