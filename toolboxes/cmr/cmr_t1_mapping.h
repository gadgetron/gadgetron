/** \file   cmr_t1_mapping.h
    \brief  Implement cardiac MR t1 mapping for 2D applications
            The input has dimension [RO E1 N S SLC]
            Temporal dimension is N
    \author Hui Xue
*/

#pragma once

#include "cmr_parametric_mapping.h"

namespace Gadgetron {

// ======================================================================================
// T1 Saturation recovery
// y = A * ( 1-exp(-ti/T1) )

template <typename T>
class CmrT1SRMapping : public CmrParametricMapping<T>
{
public:

    typedef CmrParametricMapping<T> BaseClass;
    typedef CmrT1SRMapping<T> Self;

    typedef typename BaseClass::ArrayType ArrayType;
    typedef typename BaseClass::ImageType ImageType;
    typedef typename BaseClass::ImageContinerType ImageContinerType;
    typedef typename BaseClass::VectorType VectorType;

    CmrT1SRMapping();
    virtual ~CmrT1SRMapping();

    // ======================================================================================
    /// parameter for t1 SR mapping
    // ======================================================================================


    // ======================================================================================
    // perform every steps
    // ======================================================================================

    /// provide initial guess for the mapping
    virtual void get_initial_guess(const VectorType& ti, const VectorType& yi, VectorType& guess);

    /// compute map values for every parameters in bi
    virtual void compute_map(const VectorType& ti, const VectorType& yi, const VectorType& guess, VectorType& bi, T& map_v);

    /// compute SD values for every parameters in bi
    virtual void compute_sd(const VectorType& ti, const VectorType& yi, const VectorType& bi, VectorType& sd, T& map_sd);

    /// two parameters, A, T1
    virtual size_t get_num_of_paras() const;

    // ======================================================================================
    /// parameter from BaseClass
    // ======================================================================================

    using BaseClass::fill_holes_in_maps_;
    using BaseClass::max_size_of_holes_;
    using BaseClass::hole_marking_value_;
    using BaseClass::compute_SD_maps_;
    using BaseClass::mask_for_mapping_;
    using BaseClass::ti_;
    using BaseClass::data_;
    using BaseClass::map_;
    using BaseClass::para_;
    using BaseClass::sd_map_;
    using BaseClass::sd_para_;
    using BaseClass::max_iter_;
    using BaseClass::max_fun_eval_;
    using BaseClass::thres_fun_;
    using BaseClass::max_map_value_;
    using BaseClass::min_map_value_;

    using BaseClass::verbose_;
    using BaseClass::debug_folder_;
    using BaseClass::perform_timing_;
    using BaseClass::gt_timer_local_;
    using BaseClass::gt_timer_;
    using BaseClass::gt_exporter_;
};

}
