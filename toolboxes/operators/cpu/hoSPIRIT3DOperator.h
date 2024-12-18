/** \file       hoSPIRIT3DOperator.h
    \brief      Implement SPIRIT 3D operator functinalities
    \author     Hui Xue
*/

#pragma once

#include "hoSPIRITOperator.h"

namespace Gadgetron {

template <typename T>
class hoSPIRIT3DOperator : public hoSPIRITOperator<T>
{
public:

    typedef hoSPIRITOperator<T> BaseClass;
    typedef typename BaseClass::ARRAY_TYPE ARRAY_TYPE;
    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;

    hoSPIRIT3DOperator(const std::vector<size_t>& dims);
    virtual ~hoSPIRIT3DOperator();

    using BaseClass::use_non_centered_fft_;
    //using BaseClass::performTiming_;
    //using BaseClass::debugFolder_;

protected:

    /// convert to image domain or back to kspace
    /// x : [RO E1 E2 ...]
    virtual void convert_to_image(const ARRAY_TYPE& x, ARRAY_TYPE& im);
    /// im : [RO E1 E2 srcCHA]
    virtual void convert_to_kspace(const ARRAY_TYPE& im, ARRAY_TYPE& x);

    using BaseClass::forward_kernel_;
    using BaseClass::adjoint_kernel_;
    using BaseClass::adjoint_forward_kernel_;
    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::coil_senMap_;

    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::kspace_dst_;
    using BaseClass::res_after_apply_kernel_;
    using BaseClass::res_after_apply_kernel_sum_over_;

    using BaseClass::fft_im_buffer_;
    using BaseClass::fft_kspace_buffer_;

    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::gt_exporter_;
};

}
