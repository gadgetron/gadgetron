
/** \file   mri_core_utility.h
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"
#include "mri_core_data.h"
#include "hoNDKLT.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// detect whether a readout has been sampled or not
    // --------------------------------------------------------------------------
    /// data: [RO E1 E2 CHA N S SLC]

    /// sampled: [E1 E2 N S SLC], if a readout is sampled, corresponding sampled location is 1; otherwise, 0
    template <typename T> EXPORTMRICORE hoNDArray<bool> detect_readout_sampling_status(const hoNDArray<T>& data);

    /// detect sampled region along E1
    template <typename T> EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<T>& data);

    /// detect sampled region along E2
    template <typename T> EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray<T>& data);

    /// zero padding resize for kspace and complex images
    /// if sizeE2<=1, 2D zero padding resize is performed
    template <typename T> EXPORTMRICORE void zero_pad_resize(const hoNDArray<T>& complexIm, size_t sizeRO, size_t sizeE1, size_t sizeE2, hoNDArray<T>& complexImResized);

    /// compute averaged kspace along N and S dimensions
    /// if count_sampling_freq == true, number of times where a line is sampled along N is counted, instread of an uniformed averaging
    template <typename T> EXPORTMRICORE void compute_averaged_data_N_S(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray<T>& res);

    /// compute KL coefficients from the data
    /// the KLT has the size [N S SLC]
    /// if average_N==true or average_S==true, data will first be averaged along N or S
    /// if coil_compression_thres>0 or compression_num_modesKept>0, the number of kept channels is determine; compression_num_modesKept has the priority if it is set
    /// for all N, S and SLC, the same number of channels is kept. This number is either set by compression_num_modesKept or automatically determined in the first KLT prepare call
    template <typename T> EXPORTMRICORE void compute_eigen_channel_coefficients(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, size_t N, size_t S, double coil_compression_thres, size_t compression_num_modesKept, std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT);

    /// apply eigen channel coefficients
    /// apply KLT coefficients to data for every N, S, and SLC
    template <typename T> EXPORTMRICORE void apply_eigen_channel_coefficients(const std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT, hoNDArray<T>& data);

    /// get the path of debug folder
    // environmental variable GADGETRON_DEBUG_FOLDER is used 
    EXPORTMRICORE void get_debug_folder_path(const std::string& debugFolder, std::string& debugFolderPath);
}
