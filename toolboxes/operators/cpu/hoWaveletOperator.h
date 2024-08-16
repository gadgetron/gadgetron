/** \file       hoWaveletOperator.h
    \brief      Implement wavelet operator for L1 regularization
    \author     Hui Xue
*/

#pragma once

#include "cpuOperatorExport.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "linearOperator.h"
#include "hoNDHarrWavelet.h"
#include "hoNDRedundantWavelet.h"

//#include "GadgetronTimer.h"
//#include "ImageIOAnalyze.h"

namespace Gadgetron { 

template <typename T> 
class EXPORTCPUOPERATOR hoWaveletOperator : public linearOperator< hoNDArray<T> >
{
public:

    typedef linearOperator< hoNDArray<T> > BaseClass;
    typedef hoNDArray<T> ARRAY_TYPE;

    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
    typedef typename realType<T>::Type value_type;

    hoWaveletOperator(const std::vector<size_t>& dims);
    virtual ~hoWaveletOperator();

    /// forward wavelet transform
    virtual void mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false) = 0;
    /// backward wavelet transform
    virtual void mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate = false) = 0;

    /// compute gradient of L1 wavelet operation
    virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate = false) = 0;

    /// compute cost value of L1 wavelet operation
    virtual REAL magnitude(ARRAY_TYPE* x) = 0;

    /// proximal operation for the L1 norm of wavelet
    virtual void proximity(hoNDArray<T>& wavCoeff, value_type thres) = 0;

    virtual bool unitary() const { return true; }

    // restore acquired kspace points to x
    virtual void restore_acquired_kspace(const ARRAY_TYPE& acquired, ARRAY_TYPE& y);
    virtual void restore_acquired_kspace(ARRAY_TYPE& y);

    // set the acquired kspace, unacquired points are set to be zero
    virtual void set_acquired_points(ARRAY_TYPE& kspace);

    /// select wavelet
    /// wav_name: db1, db2, db3, db4, db5
    virtual void select_wavelet(const std::string& wav_name);

    /// if true, input is in kspace domain
    bool input_in_kspace_;

    /// if input_in_kspace_==true, indicate whether null-space constraint is applied or not
    bool no_null_space_;

    // number of transformation levels
    size_t num_of_wav_levels_;

    // whether to include low frequency approximation coefficients
    bool with_approx_coeff_;

    // whether to perform proximity across channels
    bool proximity_across_cha_;

    /// if set, the coil map is used during operation
    hoNDArray<T> coil_map_;

    /// whether to perform timing
    bool performTiming_;

    /// debug folder
    std::string debugFolder_;

protected:

    /// wavelet object used
    hoNDWavelet<T>* p_active_wav_;

    /// wavelet object
    hoNDHarrWavelet<T> harr_wav_;

    /// general redundant wavelet object
    hoNDRedundantWavelet<T> redundant_wav_;

    /// for null-space type operators
    hoNDArray<T> acquired_points_;
    // acquired point indicator array, acquired points as 1, otherwise, 0
    hoNDArray<T> acquired_points_indicator_;
    // unacquired point indicator array
    hoNDArray<T> unacquired_points_indicator_;

    // helper memory
    hoNDArray<T> complexIm_;
    hoNDArray<T> complexIm_after_apply_coil_map_;
    hoNDArray<T> kspace_;

    hoNDArray<T> complexIm_fft_;
    hoNDArray<T> kspace_fft_;

    hoNDArray<T> wav_coeff_;
    hoNDArray<value_type> wav_coeff_norm_;
    hoNDArray<value_type> wav_coeff_norm_approx_;
    hoNDArray<value_type> wav_coeff_norm_mag_;

    hoNDArray<T> kspace_wav_;
    hoNDArray<T> complexIm_wav_;
    hoNDArray<value_type> complexIm_norm_;

    // clock for timing
    //Gadgetron::GadgetronTimer gt_timer1_;
    //Gadgetron::GadgetronTimer gt_timer2_;
    //Gadgetron::GadgetronTimer gt_timer3_;

    // exporter
    //Gadgetron::ImageIOAnalyze gt_exporter_;
};

}
