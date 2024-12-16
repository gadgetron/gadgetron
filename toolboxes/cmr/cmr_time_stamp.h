/** \file   cmr_time_stamp.h
    \brief  Implement functionalities to handle cardiac time stamps
    \author Hui Xue
*/

#pragma once

#include "GadgetronTimer.h"

#include "mri_core_def.h"
#include "mri_core_utility.h"
#include "hoNDBSpline.h"

namespace Gadgetron {

    /// given the incoming acquisition time stamp in ms, correct the digitizing error using a linear fit
    /// unacquired lines are marked with time_stamp <0
    /// time_stamp : [E1 N] array
    void correct_time_stamp_with_fitting(hoNDArray<float>& time_stamp, size_t startE1, size_t endE1);

    /// detect heart beat
    /// cpt_time_stamp : cardiac phase time array, [E1 N], the missing lines are marked with cpt time stamp <0
    /// ind_hb: [E1 N] array, mark which heart beat a readout line belongs; first heart beat is 0
    /// start_e1_hb, end_e1_hb: for every detected heart beat, mark its staring and ending e1 index
    /// start_n_hb, end_n_hb: for every detected heart beat, mark its staring and ending n index
    void detect_heart_beat_with_time_stamp(hoNDArray<float>& cpt_time_stamp, hoNDArray<int>& ind_hb,
                                        std::vector<size_t>& start_e1_hb, std::vector<size_t>& end_e1_hb,
                                        std::vector<size_t>& start_n_hb, std::vector<size_t>& end_n_hb );

    /// correct cpt time stamp with line fit to every heart beat
    void correct_heart_beat_time_stamp_with_fitting(hoNDArray<float>& cpt_time_stamp, hoNDArray<int>& ind_hb, size_t startE1, size_t endE1,
                                                            const std::vector<size_t>& start_e1_hb, const std::vector<size_t>& end_e1_hb,
                                                            const std::vector<size_t>& start_n_hb, const std::vector<size_t>& end_n_hb );

    /// given the filled time stamp arrays [E1 N], compute time stamp for every n
    void compute_phase_time_stamp(const hoNDArray<float>& time_stamp, const hoNDArray<float>& cpt_time_stamp, size_t startE1, size_t endE1,
        hoNDArray<float>& phs_time_stamp, hoNDArray<float>& phs_cpt_time_stamp);

    /// resample data array along last dimension, for cmr interpolation
    /// the data has the size [... N], and assumed to be sampled at [0 1 2 ... N-1]
    /// resampled res array will be computed at [0 dn 2*dn ... output_N-1], dn = (N-1)/(output_N-1)
    /// BSpline will be used for interpolation
    template <typename T> void resample_cardiac_phase_cmr_array(const hoNDArray<T>& data, size_t output_N, hoNDArray<T>& res, size_t spline_degree=5);
}
