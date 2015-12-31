/** \file       hoSPIRITOperator.h
    \brief      Implement SPIRIT operator functinalities
    \author     Hui Xue
*/

#pragma once

#include "cpuOperatorExport.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "linearOperator.h"

#include "GadgetronTimer.h"
#include "gtPlusIOAnalyze.h"

namespace Gadgetron { 

template <typename T> 
class EXPORTCPUOPERATOR hoSPIRITOperator : public linearOperator< hoNDArray<T> >
{
public:

    typedef linearOperator< hoNDArray<T> > BaseClass;
    typedef hoNDArray<T> ARRAY_TYPE;

    typedef typename BaseClass::REAL REAL;
    typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;

    hoSPIRITOperator(std::vector<size_t> *dims);
    virtual ~hoSPIRITOperator();

    /// set forward kernel, compute the adjoint and adjoint_forward kernel
    void set_forward_kernel(hoNDArray<T>& forward_kernel, bool compute_adjoint_forward_kernel=false);

    /// apply(G-I)Dc'
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

    /// indicate the operator is unitary or not
    /// unitary operator, AA' = I
    virtual bool unitary() const { return false; }

    /// convert to image domain or back to kspace
    virtual void convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im) = 0;
    virtual void convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x) = 0;

    // restore acquired kspace points to x
    virtual void restore_acquired_kspace(const hoNDArray<T>& acquired, hoNDArray<T>& y);
    virtual void restore_acquired_kspace(hoNDArray<T>& y);

    // set the acquired kspace, unacquired points are set to be zero
    virtual void set_acquired_points(hoNDArray<T>& kspace);

    // set the coil sensivity map
    virtual void set_coil_sen_map(hoNDArray<T>& senMap);

    /// if true, use the fft. not fftc
    bool use_non_centered_fft_;

    /// whether to perform timing
    bool performTiming_;

    /// debug folder
    std::string debugFolder_;

protected:

    // G-I, [... srcCHA dstCHA]
    ARRAY_TYPE forward_kernel_;
    // (G-I)', [... dstCHA srcCHA]
    ARRAY_TYPE adjoint_kernel_;
    // (G-I)'(G-I)
    ARRAY_TYPE adjoint_forward_kernel_;

    ARRAY_TYPE acquired_points_;
    // acquired point indicator array, acquired points as 1, otherwise, 0
    ARRAY_TYPE acquired_points_indicator_;
    // unacquired point indicator array
    ARRAY_TYPE unacquired_points_indicator_;

    ARRAY_TYPE coil_senMap_;

    // utility functions
    void sum_over_src_channel(const hoNDArray<T>& x, hoNDArray<T>& r);

    // helper memory
    ARRAY_TYPE kspace_;
    ARRAY_TYPE complexIm_;
    ARRAY_TYPE kspace_dst_;
    ARRAY_TYPE res_after_apply_kernel_;
    ARRAY_TYPE res_after_apply_kernel_sum_over_;

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;
};

}
