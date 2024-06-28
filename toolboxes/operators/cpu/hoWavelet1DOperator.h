/** \file       hoWavelet1DOperator.h
    \brief      Implement wavelet operator for L1 regularization
    \author     Hui Xue
*/

#pragma once

#include "hoWaveletOperator.h"

namespace Gadgetron { 

template <typename T> 
class EXPORTCPUOPERATOR hoWavelet1DBaseOperator : public hoWaveletOperator<T>
{
public:

    typedef hoWaveletOperator<T> BaseClass;
    typedef hoNDArray<T> ARRAY_TYPE;

    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
    typedef typename realType<T>::Type value_type;

    hoWavelet1DBaseOperator(const std::vector<size_t>& dims);
    virtual ~hoWavelet1DBaseOperator();

    /// perform operation 1 ||WF'x||1,    if input_in_kspace_==true and no_null_space_ is false
    /// perform operation 2 ||Wx||1,      if input_in_kspace_==false
    /// perform operation 3 ||WF'(Dc'x+D'a)||1,    if input_in_kspace_==true and and no_null_space_ is true
    /// the coil_map_ is not used here, but can be added later if needed

    /// forward wavelet transform
    /// x: [RO CHA ...]
    /// y: [RO W CHA ...], W=1+level
    /// implement Wx, or WF'x or WF'(Dc'x+D'a)
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false) = 0;
    /// backward wavelet transform
    /// x: [RO W CHA ...]
    /// y: [RO CHA ...]
    /// implement W'x, or FW'x or Dc*FW'x
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false) = 0;

    /// compute gradient of L1 wavelet operation
    /// gradient of ||Wx||1 ~= W'(Wx./((Wx)'.*(Wx)+mu)^0.5)
    /// compute graident of ||Wx||1, or ||WF'x||1 or ||WF'(Dc'x+D'a)||1
    virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate = false);

    /// compute cost value of L1 wavelet operation of ||Wx||1, or ||WF'x||1 or ||WF'(Dc'x+D'a)||1
    virtual REAL magnitude(ARRAY_TYPE* x);

    /// proximal operation for the L1 norm of wavelet, the joint sparsity across CHA is used
    virtual void proximity(hoNDArray<T>& wavCoeff, value_type thres);

    /// if set, the wavelet coefficents will be scaled by this array
    ARRAY_TYPE wav_coeff_scaling_;

    using BaseClass::input_in_kspace_;
    using BaseClass::no_null_space_;
    using BaseClass::num_of_wav_levels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::coil_map_;
    using BaseClass::performTiming_;
    using BaseClass::debugFolder_;

protected:

    // compute L1 norm of wavelet coefficients across CHA
    // waveCoeff: [RO W CHA ...], W is the wavelet coefficient dimension (e.g. for 1 level wavelet decomposition, W=2 for 1D)
    // the W=1 wavelet coefficient is the most low frequent coefficients
    void L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<value_type>& wavCoeffNorm);

    // soft-threshold or shrink the wavelet coefficients
    // the really applied threshold is mask.*thres
    void shrink_wav_coeff(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type thres, bool processApproxCoeff = false);

    // devide the wavelet coeff by norm
    void divide_wav_coeff_by_norm(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type mu, value_type p, bool processApproxCoeff = false);

    // helper functions
    virtual void forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y);
    virtual void adjoint_wav(const hoNDArray<T>& x, hoNDArray<T>& y);

    using BaseClass::p_active_wav_;

    ARRAY_TYPE adj_x_;

    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::complexIm_after_apply_coil_map_;
    using BaseClass::kspace_;

    using BaseClass::complexIm_;
    using BaseClass::wav_coeff_;
    using BaseClass::wav_coeff_norm_;
    using BaseClass::wav_coeff_norm_approx_;
    using BaseClass::wav_coeff_norm_mag_;

    using BaseClass::kspace_wav_;
    using BaseClass::complexIm_wav_;
    using BaseClass::complexIm_norm_;

    using BaseClass::complexIm_fft_;
    using BaseClass::kspace_fft_;

    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::gt_exporter_;
};

template <typename T> 
class EXPORTCPUOPERATOR hoWavelet1DOperator : public hoWavelet1DBaseOperator<T>
{
public:

    typedef hoWavelet1DBaseOperator<T> BaseClass;
    typedef hoNDArray<T> ARRAY_TYPE;

    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
    typedef float value_type;

    hoWavelet1DOperator(const std::vector<size_t>& dims);
    virtual ~hoWavelet1DOperator();

    /// forward wavelet transform
    /// x: [RO CHA ...]
    /// y: [RO W CHA ...], W=1+level
    /// implement Wx, or WF'x or WF'(Dc'x+D'a)
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);
    /// backward wavelet transform
    /// x: [RO W CHA ...]
    /// y: [RO CHA ...]
    /// implement W'x, or FW'x or Dc*FW'x
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

    using BaseClass::input_in_kspace_;
    using BaseClass::wav_coeff_scaling_;
    using BaseClass::no_null_space_;
    using BaseClass::num_of_wav_levels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::coil_map_;
    using BaseClass::performTiming_;
    using BaseClass::debugFolder_;

protected:

    virtual void convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual void convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x);

    using BaseClass::p_active_wav_;
    using BaseClass::adj_x_;
    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::complexIm_after_apply_coil_map_;
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::wav_coeff_;
    using BaseClass::wav_coeff_norm_;
    using BaseClass::wav_coeff_norm_approx_;
    using BaseClass::wav_coeff_norm_mag_;
    using BaseClass::kspace_wav_;
    using BaseClass::complexIm_wav_;
    using BaseClass::complexIm_norm_;
    using BaseClass::complexIm_fft_;
    using BaseClass::kspace_fft_;
    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::gt_exporter_;
};

template <typename T>
class EXPORTCPUOPERATOR hoWavelet1DREALOperator : public hoWavelet1DBaseOperator<T>
{
public:

    typedef hoWavelet1DBaseOperator<T> BaseClass;
    typedef hoNDArray<T> ARRAY_TYPE;

    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
    typedef float value_type;

    hoWavelet1DREALOperator(const std::vector<size_t>& dims);
    virtual ~hoWavelet1DREALOperator();

    /// forward wavelet transform
    /// x: [RO CHA ...]
    /// y: [RO W CHA ...], W=1+level
    /// implement Wx
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);
    /// backward wavelet transform
    /// x: [RO W CHA ...]
    /// y: [RO CHA ...]
    /// implement W'x
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false);

    using BaseClass::input_in_kspace_;
    using BaseClass::wav_coeff_scaling_;
    using BaseClass::no_null_space_;
    using BaseClass::num_of_wav_levels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::coil_map_;
    using BaseClass::performTiming_;
    using BaseClass::debugFolder_;

protected:

    using BaseClass::p_active_wav_;
    using BaseClass::adj_x_;
    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::complexIm_after_apply_coil_map_;
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::wav_coeff_;
    using BaseClass::wav_coeff_norm_;
    using BaseClass::wav_coeff_norm_approx_;
    using BaseClass::wav_coeff_norm_mag_;
    using BaseClass::kspace_wav_;
    using BaseClass::complexIm_wav_;
    using BaseClass::complexIm_norm_;
    using BaseClass::complexIm_fft_;
    using BaseClass::kspace_fft_;
    //using BaseClass::gt_timer1_;
    //using BaseClass::gt_timer2_;
    //using BaseClass::gt_timer3_;
    //using BaseClass::gt_exporter_;
};

}
