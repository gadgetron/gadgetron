/** \file   cmr_parametric_mapping.h
    \brief  Implement cardiac MR parametric mapping (T1/T2/T2*) for 2D applications
            The input has dimension [RO E1 N S SLC]
            Temporal dimension is N
    \author Hui Xue
*/

#pragma once

#include "GadgetronTimer.h"

#include "ImageIOAnalyze.h"
#include "mri_core_def.h"
#include "mri_core_utility.h"
#include "hoNDImageContainer2D.h"
#include "hoMRImage.h"

namespace Gadgetron {

    /// map: the 2D map for hole filling; holes is marked by value 'hole'
    /// hole: this value marks the hole
    /// is_8_connected: whethe to use 8-connection to detect holes; if false, 4-connection is used
    template <typename T>
    void perform_hole_filling(hoNDArray<T>& map, T hole=(T)(0), size_t max_size_of_holes=20, bool is_8_connected=true);
    /// hole filling on image container
    template <typename T>
    void perform_hole_filling(Gadgetron::hoNDImageContainer2D< hoMRImage<T, 2> >& maps, T hole = (T)(0), size_t max_size_of_holes = 20, bool is_8_connected = true);

    template <typename T>
    class CmrParametricMapping
    {
    public:

        typedef CmrParametricMapping<T> Self;
        typedef hoNDArray< T > ArrayType;

        // image type
        typedef Gadgetron::hoMRImage<T, 2> ImageType;
        typedef Gadgetron::hoNDImageContainer2D<ImageType> ImageContinerType;

        // vector type
        typedef std::vector<T> VectorType;

        CmrParametricMapping();
        virtual ~CmrParametricMapping();

        // ======================================================================================
        // perform the mapping
        // ======================================================================================
        virtual void perform_parametric_mapping();

        // ======================================================================================
        /// parameter for overall workflow
        // ======================================================================================

        // whether true, hole filling will be performed on maps
        bool fill_holes_in_maps_;
        /// maximal size for holes
        size_t max_size_of_holes_;
        /// marking value for holes in maps
        T hole_marking_value_;

        /// whether to compute SD maps
        bool compute_SD_maps_;

        /// mask for mapping, pixels used for mapping is marked as >0
        /// if empty, every pixel is inputted for mapping
        hoNDArray<T> mask_for_mapping_;

        // ======================================================================================
        /// parameter for mapping
        // ======================================================================================

        /// parametric times
        /// inversion/saturation/TE times for every image along N
        std::vector<T> ti_;

        /// input array, [RO E1 N S SLC]
        hoNDArray<T> data_;

        /// map array, [RO E1 S SLC]
        hoNDArray<T> map_;

        /// parameters array, [RO E1 NUM S SLC]
        /// there are NUM parameters (including the map)
        hoNDArray<T> para_;

        /// sd map array, [RO E1 S SLC]
        hoNDArray<T> sd_map_;

        /// sd array for other parameters, [RO E1 NUM S SLC]
        hoNDArray<T> sd_para_;

        /// common parameters for optimization
        /// maximal number of iterations
        size_t max_iter_;
        /// maximal number of function evaluation
        size_t max_fun_eval_;
        /// threshold for minimal function value change
        T thres_fun_;

        /// maximal valid value of map
        T max_map_value_;
        T min_map_value_;

        // ======================================================================================
        /// parameter for debugging
        // ======================================================================================

        // verbose mode for more output messages
        bool verbose_;
        // debug folder
        std::string debug_folder_;
        // whether to perform timing
        bool perform_timing_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // ======================================================================================
        // perform every steps
        // ======================================================================================

        /// provide initial guess for the mapping
        virtual void get_initial_guess(const VectorType& ti, const VectorType& yi, VectorType& guess);

        /// compute map values for every parameters in bi
        virtual void compute_map(const VectorType& ti, const VectorType& yi, const VectorType& guess, VectorType& bi, T& map_v);

        /// compute SD values for every parameters in bi
        virtual void compute_sd(const VectorType& ti, const VectorType& yi, const VectorType& bi, VectorType& sd, T& map_sd);

        /// compute SD values from gradient vector
        virtual void compute_sd_impl(const VectorType& ti, const VectorType& yi, const VectorType& bi, const VectorType& res, const hoNDArray<T>& grad, VectorType& sd);

        /// return number of parameters, including the map itself
        virtual size_t get_num_of_paras() const;
    };
}
