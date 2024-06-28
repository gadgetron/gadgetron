/** \file       hoWavelet3DOperator.h
    \brief      Implement wavelet operator for L1 regularization
    \author     Hui Xue
*/

#pragma once

#include "hoWavelet2DTOperator.h"

namespace Gadgetron { 

template <typename T> 
class EXPORTCPUOPERATOR hoWavelet3DOperator : public hoWavelet2DTOperator<T>
{
public:

    typedef hoWavelet2DTOperator<T> BaseClass;
    typedef hoNDArray<T> ARRAY_TYPE;

    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
    typedef typename realType<T>::Type value_type;

    hoWavelet3DOperator(const std::vector<size_t>& dims);
    virtual ~hoWavelet3DOperator();

    /// if no_null_space_ == true
    /// perform operation ||WS'F'x||1,    if input_in_kspace_==true  and coil_map_ is not empty
    /// perform operation ||WS'x||1,      if input_in_kspace_==false and coil_map_ is not empty
    /// perform operation ||Wx||1,        if input_in_kspace_==false and coil_map_ is empty
    /// perform operation ||WF'x||1,      if input_in_kspace_==true  and coil_map_ is empty

    /// if no_null_space_ == true
    /// perform operation ||WS'F'(Dc'x+D'a)||1,    if coil_map_ is not empty
    /// perform operation ||WF'(Dc'x+D'a)||1,      if coil_map_ is empty

    /// forward wavelet transform operation 1/2/3/4
    /// x: [RO E1 E2 CHA]
    /// y: [RO E1 E2 W CHA], W=1+7*level
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);
    /// backward wavelet transform for operation 1/2/3/4
    /// x: [RO E1 E2 W CHA]
    /// y: [RO E1 E2 CHA]
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

    using BaseClass::input_in_kspace_;
    using BaseClass::no_null_space_;
    using BaseClass::scale_factor_first_dimension_;
    using BaseClass::scale_factor_second_dimension_;
    using BaseClass::scale_factor_third_dimension_;
    using BaseClass::change_coeffcients_third_dimension_boundary_;
    using BaseClass::num_of_wav_levels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::coil_map_;

    using BaseClass::performTiming_;
    using BaseClass::debugFolder_;

protected:

    // helper functions
    virtual void convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual void convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x);

    virtual void forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y);
    virtual void adjoint_wav(const hoNDArray<T>& x, hoNDArray<T>& y);

    using BaseClass::mask_;
    using BaseClass::forward_buf_;
    using BaseClass::adjoint_buf_;

    using BaseClass::p_active_wav_;

    using BaseClass::complexIm_;
    using BaseClass::complexIm_after_apply_coil_map_;
    using BaseClass::kspace_;

    using BaseClass::complexIm_fft_;
    using BaseClass::kspace_fft_;

    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;

    using BaseClass::wav_coeff_;
    using BaseClass::wav_coeff_norm_;
    using BaseClass::wav_coeff_norm_approx_;
    using BaseClass::wav_coeff_norm_mag_;

    using BaseClass::kspace_wav_;
    using BaseClass::complexIm_wav_;
    using BaseClass::complexIm_norm_;

    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::gt_exporter_;
};

}
