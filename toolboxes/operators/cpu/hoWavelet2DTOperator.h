/** \file       hoWavelet2DTOperator.h
    \brief      Implement wavelet operator for L1 regularization
    \author     Hui Xue
*/

#pragma once

#include "hoWaveletOperator.h"
#include "hoMotionCompensation2DTOperator.h"

namespace Gadgetron { 

template <typename T> 
class EXPORTCPUOPERATOR hoWavelet2DTOperator : public hoWaveletOperator<T>
{
public:

    typedef hoWaveletOperator<T> BaseClass;
    typedef hoNDArray<T> ARRAY_TYPE;

    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
    typedef typename realType<T>::Type value_type;

    hoWavelet2DTOperator(const std::vector<size_t>& dims);
    virtual ~hoWavelet2DTOperator();

    /// if motion componsenation is not used
        /// if no_null_space_ == true
        /// perform operation ||WS'F'x||1,    if input_in_kspace_==true  and coil_map_ is not empty
        /// perform operation ||WS'x||1,      if input_in_kspace_==false and coil_map_ is not empty
        /// perform operation ||Wx||1,        if input_in_kspace_==false and coil_map_ is empty
        /// perform operation ||WF'x||1,      if input_in_kspace_==true  and coil_map_ is empty

        /// if no_null_space_ == false
        /// perform operation ||WS'F'(Dc'x+D'a)||1,    if coil_map_ is not empty
        /// perform operation ||WF'(Dc'x+D'a)||1,      if coil_map_ is empty

    /// if motion componsenation is used
        /// if no_null_space_ == true
        /// perform operation ||WMS'F'x||1,    if input_in_kspace_==true  and coil_map_ is not empty
        /// perform operation ||WMS'x||1,      if input_in_kspace_==false and coil_map_ is not empty
        /// perform operation ||WMx||1,        if input_in_kspace_==false and coil_map_ is empty
        /// perform operation ||WMF'x||1,      if input_in_kspace_==true  and coil_map_ is empty

        /// if no_null_space_ == false
        /// perform operation ||WMS'F'(Dc'x+D'a)||1,    if coil_map_ is not empty
        /// perform operation ||WMF'(Dc'x+D'a)||1,      if coil_map_ is empty

    /// forward wavelet transform operation
    /// x: [RO E1 CHA N]
    /// y: [RO E1 N W CHA] or [RO E1 N W 1], W=1+7*level
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);
    /// backward wavelet transform for operation
    /// x: [RO E1 N W CHA] or [RO E1 N W 1]
    /// y: [RO E1 CHA N]
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

    /// compute gradient of L1 wavelet operation, depending on the value of input_in_kspace_ and coil_map_
    virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate = false);

    /// compute cost value of L1 wavelet operation, depending on the value of input_in_kspace_ and coil_map_
    virtual REAL magnitude(ARRAY_TYPE* x);

    /// proximal operation for the L1 norm of wavelet, the joint sparsity across CHA is used
    virtual void proximity(hoNDArray<T>& wavCoeff, value_type thres);

    // because the spatial resolution of images are often different in through-plane dimension than the other two dimensions
    // sometime it is good to take this into account, so the regularization effects are more isotropic
    // Here only simple scaling factors are used
    // More generally, a weighting matrix can be concatenated with wavelet coefficients to enhance or suppress regularization effects as needed
    // the regularization term can become ||Q*W*F'*(Dc'x+D'y)||1, Q is the general weighting matrix
    // in the next version, we shall extend this class with more geneal weighting strategy
    // these scaling factors are only used for high frequency wavelet coefficients!!!
    value_type scale_factor_first_dimension_;
    value_type scale_factor_second_dimension_;
    value_type scale_factor_third_dimension_;

    // in some cases, the boundary high frequency coefficients of the 3rd dimension should not be changed, e.g. the dynamic contrast experiments
    bool change_coeffcients_third_dimension_boundary_;

    /// MOCO operator
    hoMotionCompensation2DTOperator<T, double> mocoer_;

    using BaseClass::input_in_kspace_;
    using BaseClass::no_null_space_;
    using BaseClass::num_of_wav_levels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::proximity_across_cha_;
    using BaseClass::coil_map_;

    using BaseClass::performTiming_;
    using BaseClass::debugFolder_;

protected:

    // compute L1 norm of wavelet coefficients across CHA
    // waveCoeff: [RO E1 N W CHA ...], W is the wavelet coefficient dimension
    // the W=1 wavelet coefficient is the most low frequent coefficients
    void L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<value_type>& wavCoeffNorm);

    // soft-threshold or shrink the wavelet coefficients
    // the really applied threshold is mask.*thres
    void shrink_wav_coeff(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, const hoNDArray<T>& mask, bool processApproxCoeff = false);

    // devide the wavelet coeff by norm
    void divide_wav_coeff_by_norm(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type mu, value_type p, bool processApproxCoeff = false);

    // helper functions
    virtual void convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual void convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x);

    virtual void forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y);
    virtual void adjoint_wav(const hoNDArray<T>& x, hoNDArray<T>& y);

    // apply scaling along dimensions for high frequency wavelet coefficients
    void apply_scale_first_dimension(hoNDArray<T>& wavCoeff, value_type& scaleFactor);
    void apply_scale_second_dimension(hoNDArray<T>& wavCoeff, value_type& scaleFactor);
    void apply_scale_third_dimension(hoNDArray<T>& wavCoeff, value_type& scaleFactor);

    hoNDArray<T> mask_;
    hoNDArray<T> forward_buf_;
    hoNDArray<T> adjoint_buf_;

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

    hoNDArray<T> moco_image_;
    hoNDArray<T> adj_moco_image_;

/*    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::gt_exporter_*/;
};

}
