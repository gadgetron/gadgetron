/** \file   cmr_time_stamp.h
    \brief  Implement functionalities to handle cardiac time stamps
    \author Hui Xue
*/

#pragma once

#include "cmr_export.h"

#include "GadgetronTimer.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

namespace Gadgetron { 

    /// given the incoming acquisition time stamp in ms, correct the digitizing error using a linear fit
    /// unacquired lines are marked with time_stamp <0
    /// time_stamp : [E1 N] array
    EXPORTCMR void correct_time_stamp_with_fitting(hoNDArray<float>& time_stamp);

    /// detect heart beat
    /// cpt_time_stamp : cardiac phase time array, [E1 N], the missing lines are marked with cpt time stamp <0
    /// ind_hb: [E1 N] array, mark which heart beat a readout line belongs; first heart beat is 0
    /// start_e1_hb, end_e1_hb: for every detected heart beat, mark its staring and ending e1 index
    /// start_n_hb, end_n_hb: for every detected heart beat, mark its staring and ending n index
    EXPORTCMR void detect_heart_beat_with_time_stamp(hoNDArray<float>& cpt_time_stamp, hoNDArray<int>& ind_hb, 
                                        std::vector<size_t>& start_e1_hb, std::vector<size_t>& end_e1_hb, 
                                        std::vector<size_t>& start_n_hb, std::vector<size_t>& end_n_hb );

    /// correct cpt time stamp with line fit to every heart beat
    EXPORTCMR void correct_heart_beat_time_stamp_with_fitting(hoNDArray<float>& cpt_time_stamp, hoNDArray<int>& ind_hb, 
                                                            const std::vector<size_t>& start_e1_hb, const std::vector<size_t>& end_e1_hb, 
                                                            const std::vector<size_t>& start_n_hb, const std::vector<size_t>& end_n_hb );

    /// given the filled time stamp arrays [E1 N], compute time stamp for every n
    EXPORTCMR void compute_phase_time_stamp(const hoNDArray<float>& time_stamp, const hoNDArray<float>& cpt_time_stamp,
        hoNDArray<float>& phs_time_stamp, hoNDArray<float>& phs_cpt_time_stamp);
}
