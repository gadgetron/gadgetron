/** \file   cmr_motion_correction.h
    \brief  Implement functionalities to perform some cardiac MR motion correction
    \author Hui Xue
*/

#pragma once

#include "cmr_export.h"

#include <boost/range/adaptor/strided.hpp>
#include <range/v3/action.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#ifdef USE_OMP
    #include <omp.h>
#endif
#ifdef max
    #undef max
#endif
#ifdef min
    #undef min
#endif

#include <algorithm>
#include "hoMatrix.h"

#include "GadgetronTimer.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

#include "hoImageRegContainer2DRegistration.h"

namespace Gadgetron { 

    struct cmr_moco_ave_compObj
    {
        cmr_moco_ave_compObj() {}
        ~cmr_moco_ave_compObj() {}

        bool operator()(const std::pair<double, size_t>& m1, const std::pair<double, size_t>& m2) const
        {
            return !(m1.first >= m2.first);
        }
    };

    /// given the series of images along row r, compute a key frame for moco, using SSD strategy
    /// input: [RO E1 N]
    template <typename T> EXPORTCMR void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<T>& input, size_t& key_frame);

    /// given the series of images along row r, compute a key frame for moco, using SSD strategy
    /// input: [RO E1 N]
    /// moco_quality is marked the quality of motion correction
    /// first element is the quality score (smaller value means better quality)
    /// second element is the frame index
    template <typename T> EXPORTCMR void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<T>& input, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality);

    /// compute ssd based moco quality, given the found key_frame
    template <typename T> EXPORTCMR void compute_SSD_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame, std::vector< std::pair<double, size_t> >& moco_quality);

    /// compute jacobian parameters from deformation fields
    template <typename T> EXPORTCMR void compute_deformation_jacobian(const Gadgetron::hoNDArray<T>& dx, const Gadgetron::hoNDArray<T>& dy, std::vector<T>& mean_deform, std::vector<T>& max_deform, std::vector<T>& mean_log_jac, std::vector<T>& max_log_jac);

    /// find the key frame from the deformation field
    /// the minimal total deformation is used
    template <typename T> EXPORTCMR void find_key_frame_deformation_2DT(const Gadgetron::hoNDArray<T>& dx, const Gadgetron::hoNDArray<T>& dy, size_t& key_frame, std::vector<
        std::pair<double, size_t> >& moco_quality);

    /// perform moco average
    template <typename T> EXPORTCMR void perform_averaging(const Gadgetron::hoNDArray<T>& input, size_t& key_frame, const std::vector< std::pair<double, size_t> >& moco_quality, double percentage_kept_for_averaging, Gadgetron::hoNDArray<T>& ave);

    /// perform 2D moco across a series
    /// input: [RO E1 N]
    /// the 2D image series are motion corrected and deformation fields are stored in reg.deformation_field_
    /// reg_strength : regularization strength in the unit of pixel
    template <typename T> EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double >& reg);

    template <typename T> EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame,
        T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double >& reg);

    /// perform 2D moco on input image container
    template <typename T> EXPORTCMR void perform_moco_fixed_key_frame_2DT(Gadgetron::hoNDImageContainer2D< hoNDImage<T, 2> >& input, const std::vector<unsigned int>& key_frame,
        T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg);

    /// series a is mocoed to series b
    /// a, b : [RO E1 N]
    template <typename T> EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& a, const Gadgetron::hoNDArray<T>& b, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg);

    template <typename T> EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& a, const Gadgetron::hoNDArray<T>& b,
        T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg);

    /// apply the deformation field
    /// input: [RO E1 N]
    /// dx, dy: [RO E1 N], deformation fields along x (first dimension) and y (second dimension)
    /// output: warpped image array
    template <typename T> EXPORTCMR void apply_deformation_field(const Gadgetron::hoNDArray<T>& input, const Gadgetron::hoNDArray<double>& dx, const Gadgetron::hoNDArray<double>& dy, Gadgetron::hoNDArray<T>& output, Gadgetron::GT_BOUNDARY_CONDITION bh = GT_BOUNDARY_CONDITION_BORDERVALUE);

    /// concatenate motion fields
    /// dx, dy: [RO E1 N], starting from key_frame, store the deformation fields from direct neighor to current frame
    /// dx_out, dy_out: [RO E1 N], store deformation fields from key_frame to current frame
    template <typename T> EXPORTCMR void concatenate_deform_fields_2DT(const hoNDArray<T>& dx, const hoNDArray<T>& dy, size_t key_frame, hoNDArray<T>& dx_out, hoNDArray<T>& dy_out);

    /// find key frame by cross registration between two series (target and source)
    /// the goal is to find key frame for target series
    /// first, every image in source series is registered to target series
    /// the resulting deformation fields are analyzed for mean/max deformation
    /// the frame of target sereis with total deformation is picked as the key frame
    /// target: [RO E1 N]; source: [RO E1 M]
    /// dx, dy: [RO E1 N M] deformation fields
    template <typename T> EXPORTCMR void find_key_frame_use_deformation_cross_series(const Gadgetron::hoNDArray<T>& target, const Gadgetron::hoNDArray<T>& source, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality, hoNDArray<double>& dx, hoNDArray<double>& dy);
}
