/**
 * \file SDC.h
 * \brief Sampling density compensation.
 */

#pragma once

#include "SDC_kernel.h"

#include "vector_td.h"

#include <algorithm>

namespace Gadgetron
{
    /**
     * \brief Convolution mode.
     */
    enum class SDC_conv_mode
    {
        C2NC,   /**< Cartesian to non-Cartesian. */
        NC2C    /**< Non-Cartesian to Cartesian. */
    };

    /**
     * \brief Implementation of sampling density compensation.
     * 
     * This class is intended for internal use only.
     * 
     * \tparam ARRAY Array type.
     * \tparam REAL Floating-point type.
     * \tparam D Number of dimensions.
     */
    template<template<class> class ARRAY, class REAL, unsigned int D>
    class SDC_impl
    {
    public:

        SDC_impl(const vector_td<size_t, D>& matrix_size, REAL os_factor, size_t num_iterations);

        virtual ~SDC_impl() = default;

        virtual boost::shared_ptr<ARRAY<REAL>> compute(const ARRAY<vector_td<REAL, D>>& traj);

        virtual boost::shared_ptr<ARRAY<REAL>> compute(const ARRAY<vector_td<REAL, D>>& traj, const ARRAY<REAL>& initial_dcw);

    protected:

        virtual void preprocess(const ARRAY<vector_td<REAL, D>>& traj);

        virtual void convolve(const ARRAY<REAL>& in, ARRAY<REAL>& out, SDC_conv_mode mode) = 0;

        virtual void update(const ARRAY<REAL>& in, ARRAY<REAL>& out) = 0;
        
        vector_td<size_t, D> matrix_size_;

        vector_td<size_t, D> grid_size_;

        REAL os_factor_;

        size_t num_iterations_;

        size_t num_samples_;

        size_t num_frames_;

        SDC_kernel<ARRAY, REAL, D> kernel_;
        
    };

    namespace SDC_internal
    {
        template<class T>
        struct safe_divides
        {
            __host__ __device__ T operator()(const T& x, const T& y) const
            {
                return y == T(0) ? T(0) : x / y;
            }
        };   
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    SDC_impl<ARRAY, REAL, D>::SDC_impl(const vector_td<size_t, D>& matrix_size, REAL os_factor, size_t num_iterations)
      : matrix_size_(matrix_size)
      , grid_size_(vector_td<size_t, D>(vector_td<REAL, D>(matrix_size) * os_factor))
      , os_factor_(os_factor)
      , num_iterations_(num_iterations)
      , num_samples_(0)
      , num_frames_(0)
    {
        
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    boost::shared_ptr<ARRAY<REAL>> SDC_impl<ARRAY, REAL, D>::compute(const ARRAY<vector_td<REAL, D>>& traj)
    {
        ARRAY<REAL> dcw(*traj.get_dimensions());
        fill(&dcw, (REAL)1);
        return compute(traj, dcw);
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    boost::shared_ptr<ARRAY<REAL>> SDC_impl<ARRAY, REAL, D>::compute(const ARRAY<vector_td<REAL, D>>& traj, const ARRAY<REAL>& initial_dcw)
    {
        assert(*traj.get_dimensions() == *initial_dcw.get_dimensions());

        // preprocess
        preprocess(traj);

        // working arrays
        ARRAY<REAL> dcw(initial_dcw);
        ARRAY<REAL> grid(to_std_vector(grid_size_));
        ARRAY<REAL> tmp(*dcw.get_dimensions());

        // iteration loop
        for (size_t i = 0; i < num_iterations_; i++)
        {
            // to intermediate grid
            convolve(dcw, grid, SDC_conv_mode::NC2C);

            // to original trajectory
            convolve(grid, tmp, SDC_conv_mode::C2NC);

            // update weights
            update(tmp, dcw);
        }

        return boost::make_shared<ARRAY<REAL>>(dcw);
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    void SDC_impl<ARRAY, REAL, D>::preprocess(const ARRAY<vector_td<REAL,D>>& traj)
    {
        num_samples_ = traj.get_size(0);
        num_frames_ = traj.get_number_of_elements() / num_samples_;
    }

}   // namespace Gadgetron