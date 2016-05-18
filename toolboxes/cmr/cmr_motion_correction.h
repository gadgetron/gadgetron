/** \file   cmr_motion_correction.h
    \brief  Implement functionalities to perform some cardiac MR motion correction
    \author Hui Xue
*/

#pragma once

#include "cmr_export.h"

#ifdef min
    #undef min
#endif // min

#include <algorithm>
#include "hoMatrix.h"

#include "GadgetronTimer.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

#include "hoImageRegContainer2DRegistration.h"

namespace Gadgetron { 

    /// given the series of images along row r, compute a key frame for moco, using SSD strategy
    /// input: [RO E1 N]
    template <typename T> EXPORTCMR void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<T>& input, size_t& key_frame);

    /// perform 2D moco across a series
    /// input: [RO E1 N]
    /// the 2D image series are motion corrected and deformation fields are stored in reg.deformation_field_
    /// reg_strength : regularization strength in the unit of pixel
    template <typename T> EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame, 
        T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg);

    /// series a is mocoed to series b
    /// a, b : [RO E1 N]
    template <typename T> EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& a, const Gadgetron::hoNDArray<T>& b, 
        size_t key_frame, T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg);
}
