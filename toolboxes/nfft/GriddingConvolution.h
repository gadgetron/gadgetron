#pragma once

#include <memory>

#include "complext.h"
#include "vector_td.h"

#include "ConvolutionKernel.h"

namespace Gadgetron
{
    /**
     * \brief Gridding convolution factory.
     *
     * This factory class is meant to be specialized for each array type.
     * 
     * \tparam ARRAY Array type.
     * \tparam T Value type. Can be real or complex.
     * \tparam D Number of dimensions.
     * \tparam K Convolution kernel type.
     */
    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    struct GriddingConvolution
    {
        
    };


    /**
     * \brief Direction of the gridding convolution.
     */
    enum class GriddingConvolutionMode
    {
        C2NC,   /**< Cartesian to non-Cartesian. */
        NC2C    /**< Non-Cartesian to Cartesian. */
    };


    /**
     * \brief Preprocessing mode.
     */
    enum class GriddingConvolutionPrepMode
    {
        C2NC,   /**< Prepare for C2NC mode only. */
        NC2C,   /**< Prepare for NC2C mode only. */
        ALL     /**< Prepare for both modes. */
    };


    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    class GriddingConvolutionBase
    {
    public:

        typedef realType_t<T> REAL;

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Size of the grid.
         * \param matrix_size_os Size of the oversampled grid.
         * \param kernel Convolution kernel.
         */
        GriddingConvolutionBase(const vector_td<size_t, D>& matrix_size,
                                const vector_td<size_t, D>& matrix_size_os,
                                const K<REAL, D>& kernel);

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Size of the grid.
         * \param os_factor Oversampling factor.
         * \param kernel Convolution kernel.
         */
        GriddingConvolutionBase(const vector_td<size_t, D>& matrix_size,
                                REAL os_factor,
                                const K<REAL, D>& kernel);

        virtual ~GriddingConvolutionBase() = default;
        
        /**
         * \brief Prepare gridding convolution.
         * 
         * \param trajectory Trajectory, normalized to [-0.5, 0.5].
         * \param prep_mode Preparation mode.
         */
        virtual void preprocess(const ARRAY<vector_td<REAL, D>> &trajectory,
                                GriddingConvolutionPrepMode prep_mode =
                                    GriddingConvolutionPrepMode::ALL);
        
        /**
         * \brief Perform convolution.
         * 
         * \param[in] input Input array.
         * \param[out] output Output array.
         * \param[in] mode GriddingConvolutionMode.
         * \param[in] accumulate If true, accumulate result to output array.
         *                       If false, overwrite output array.
         */
        virtual void compute(const ARRAY<T>& input,
                             ARRAY<T>& output,
                             GriddingConvolutionMode mode,
                             bool accumulate = false);

        /**
         * \brief Get matrix size.
         */
        vector_td<size_t, D> get_matrix_size() const;

        /**
         * \brief Get matrix size with oversampling.
         */
        vector_td<size_t, D> get_matrix_size_os() const;

        /**
         * \brief Get a constant reference to convolution kernel.
         */
        const K<REAL, D>& get_kernel() const;

        /**
         * \brief Get number of non-Cartesian samples.
         */
        size_t get_num_samples() const;

        /**
         * \brief Get number of frames.
         */
        size_t get_num_frames() const;

        
    protected:

        /**
         * \brief Compute gridding convolution (Cartesian to non-Cartesian).
         * 
         * \param[in] image Image.
         * \param[out] samples Non-Cartesian samples.
         * \param[in] accumulate If true, accumulate the result to output array.
         *                       If false, overwrite output array.
         */
        virtual void compute_C2NC(const ARRAY<T>& image,
                                  ARRAY<T>& samples,
                                  bool accumulate) = 0;

        /**
         * \brief Compute gridding convolution (non-Cartesian to Cartesian).
         * 
         * \param[in] samples Non-Cartesian samples.
         * \param[out] image Image.
         * \param[in] accumulate If true, accumulate the result to output array.
         *                       If false, overwrite output array.
         */
        virtual void compute_NC2C(const ARRAY<T>& samples,
                                  ARRAY<T>& image,
                                  bool accumulate) = 0;

        vector_td<size_t, D> matrix_size_;

        vector_td<size_t, D> matrix_size_os_;

        K<REAL, D> kernel_;
      
        size_t num_samples_;

        size_t num_frames_;
    };

} // namespace Gadgetron

#include "GriddingConvolution.hpp"
