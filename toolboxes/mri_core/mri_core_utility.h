
/** \file   mri_core_utility.h
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"
#include "mri_core_data.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// detect whether a readout has been sampled or not
    // --------------------------------------------------------------------------
    /// data: [RO E1 E2 CHA N S SLC]

    /// sampled: [E1 E2 N S SLC], if a readout is sampled, corresponding sampled location is 1; otherwise, 0
    template <typename T> EXPORTMRICORE void detect_readout_sampling_status(const hoNDArray<T>& data, hoNDArray<float>& sampled);

    /// detect sampled region along E1
    template <typename T> EXPORTMRICORE void detect_sampled_region_E1(const hoNDArray<T>& data, size_t& start_E1, size_t& end_E1);

    /// detect sampled region along E2
    template <typename T> EXPORTMRICORE void detect_sampled_region_E2(const hoNDArray<T>& data, size_t& start_E2, size_t& end_E2);

    /// set up the kspace filter for ref used for coil map estimation
    template <typename T> EXPORTMRICORE void generate_ref_filter_for_coil_map(const hoNDArray<T>& ref, const SamplingLimit& lim_RO, const SamplingLimit& lim_E1, const SamplingLimit& lim_E2, hoNDArray<T>& filter_RO, hoNDArray<T>& filter_E1, hoNDArray<T>& filter_E2);

}
