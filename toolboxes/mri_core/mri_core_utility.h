
/** \file   mri_core_utility.h
    \brief  Implementation useful utility functionalities for 2D and 3D MRI reconstruction
            The data structures defined in mri_core_data.h are used
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"
#include "mri_core_data.h"
#include "mri_core_kspace_filter.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// detect whether a readout has been sampled or not
    // --------------------------------------------------------------------------
    /// data: [RO E1 E2 CHA N S SLC]

    /// sampled: [E1 E2 N S SLC], if a readout is sampled, corresponding sampled location is 1; otherwise, 0
    template <typename T> EXPORTMRICORE void detect_readout_sampling_status(const hoNDArray<T>& data, hoNDArray<float>& sampled );

    /// average kspace across N
    template <typename T> EXPORTMRICORE void average_kspace_across_N(const hoNDArray<T>& data, hoNDArray<T>& averaged);

    /// average kspace across S
    template <typename T> EXPORTMRICORE void average_kspace_across_S(const hoNDArray<T>& data, hoNDArray<T>& averaged);

    /// detect sampled region along E1
    template <typename T> EXPORTMRICORE void detect_sampled_region_E1(const hoNDArray<T>& data, size_t& start_E1, size_t& end_E1);

    /// detect sampled region along E2
    template <typename T> EXPORTMRICORE void detect_sampled_region_E2(const hoNDArray<T>& data, size_t& start_E2, size_t& end_E2);

    /// detect sampled region along E2

    /// set up the kspace filter for ref used for coil map estimation
    template <typename T> EXPORTMRICORE void generate_ref_filter_for_coil_map(const hoNDArray<T>& ref, const SamplingLimit& lim_RO, const SamplingLimit& lim_E1, const SamplingLimit& lim_E2, hoNDArray<T>& filter_RO, hoNDArray<T>& filter_E1, hoNDArray<T>& filter_E2);

    /// zero padding resize for kspace and complex images
    /// if sizeE2<=1, 2D zero padding resize is performed
    template <typename T> EXPORTMRICORE void zero_pad_resize(const hoNDArray<T>& complexIm, size_t sizeRO, size_t sizeE1, size_t sizeE2, hoNDArray<T>& complexImResized);

    /// get the path of debug folder
    // environmental variable GADGETRON_DEBUG_FOLDER is used 
    EXPORTMRICORE void get_debug_folder_path(const std::string& debugFolder, std::string& debugFolderPath);
}
