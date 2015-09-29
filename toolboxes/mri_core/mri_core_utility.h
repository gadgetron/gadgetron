
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
    template <typename T> EXPORTMRICORE std::tuple<hoNDArray<bool> > detect_readout_sampling_status(const hoNDArray<T>& data);

    /// detect sampled region along E1
    template <typename T> EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<T>& data);

    /// detect sampled region along E2
    template <typename T> EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray<T>& data);
}
