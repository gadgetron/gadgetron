/** \file       hoSPIRIT2DTOperator.h
    \brief      Implement SPIRIT operator functinalities for 2D+T cases
    \author     Hui Xue
*/

#pragma once

#include "hoSPIRITOperator.h"

namespace Gadgetron { 

template <typename T> 
class EXPORTCPUOPERATOR hoSPIRIT2DTOperator : public hoSPIRITOperator<T>
{
public:

    typedef hoSPIRITOperator<T> BaseClass;
    typedef typename BaseClass::ARRAY_TYPE ARRAY_TYPE;
    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;

    hoSPIRIT2DTOperator(const std::vector<size_t>& dims);
    virtual ~hoSPIRIT2DTOperator();

    /// set forward kernel, compute the adjoint and adjoint_forward kernel
    /// forward_kernel : [RO E1 srcCHA dstCHA Nor1]
    /// the number of kernels can be 1 or equal to the number of 2D kspaces
    void set_forward_kernel(ARRAY_TYPE& forward_kernel, bool compute_adjoint_forward_kernel=false);

    /// apply(G-I)Dc'
    /// x: [RO E1 srcCHA N]
    /// y: [RO E1 dstCHA N]
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);
    /// apply Dc(G-I)'
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

    /// compute right hand side
    /// b = -(G-I)D'x
    virtual void compute_righ_hand_side(const ARRAY_TYPE& x, ARRAY_TYPE& b);

    /// compute gradient of ||(G-I)(Dc'x+D'y)||2
    virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate = false);

    /// compute cost value of L2 norm ||(G-I)(Dc'x+D'y)||2
    virtual REAL magnitude(ARRAY_TYPE* x);

    using BaseClass::use_non_centered_fft_;
    using BaseClass::no_null_space_;
    //using BaseClass::performTiming_;
    //using BaseClass::debugFolder_;

protected:

    /// convert to image domain or back to kspace
    virtual void convert_to_image(const ARRAY_TYPE& x, ARRAY_TYPE& im);
    virtual void convert_to_kspace(const ARRAY_TYPE& im, ARRAY_TYPE& x);

    using BaseClass::forward_kernel_;
    using BaseClass::adjoint_kernel_;
    using BaseClass::adjoint_forward_kernel_;
    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::coil_senMap_;

    using BaseClass::kspace_;
    using BaseClass::kspace_dst_;
    using BaseClass::complexIm_;
    ARRAY_TYPE complexIm_dst_;
    using BaseClass::res_after_apply_kernel_;
    ARRAY_TYPE res_after_apply_kernel_dst_;
    using BaseClass::res_after_apply_kernel_sum_over_;
    ARRAY_TYPE res_after_apply_kernel_sum_over_dst_;

    using BaseClass::fft_im_buffer_;
    using BaseClass::fft_kspace_buffer_;

    ARRAY_TYPE fft_im_buffer_dst_;
    ARRAY_TYPE fft_kspace_buffer_dst_;

    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::gt_exporter_;

    void apply_forward_kernel(ARRAY_TYPE& x);
    void apply_adjoint_kernel(ARRAY_TYPE& x);
    void apply_adjoint_forward_kernel(ARRAY_TYPE& x);
};

}
