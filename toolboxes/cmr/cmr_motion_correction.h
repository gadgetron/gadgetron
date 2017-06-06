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
        T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg);
    /// perform 2D moco on input image container
    template <typename T> EXPORTCMR void perform_moco_fixed_key_frame_2DT(Gadgetron::hoNDImageContainer2D< hoNDImage<T, 2> >& input, const std::vector<unsigned int>& key_frame, 
        T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg);

    /// series a is mocoed to series b
    /// a, b : [RO E1 N]
    template <typename T> EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& a, const Gadgetron::hoNDArray<T>& b, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg);

    template <typename T> EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& a, const Gadgetron::hoNDArray<T>& b, 
        T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg);

    /// apply the deformation field
    /// input: [RO E1 N]
    /// dx, dy: [RO E1 N], deformation fields along x (first dimension) and y (second dimension)
    /// output: warpped image array
    template <typename T> EXPORTCMR void apply_deformation_field(const Gadgetron::hoNDArray<T>& input, const Gadgetron::hoNDArray<float>& dx, const Gadgetron::hoNDArray<float>& dy, Gadgetron::hoNDArray<T>& output, Gadgetron::GT_BOUNDARY_CONDITION bh=GT_BOUNDARY_CONDITION_BORDERVALUE);
}
