/**
 * \file cuSDC.h
 * \brief Sampling density compensation (CUDA specialization).
 */

#pragma once

#include "SDC.h"
#include "mri_sdc_export.h"

#include "cuNDArray.h"
#include "vector_td.h"

#include <thrust/device_vector.h>


namespace Gadgetron
{
    /**
     * \brief Implementation of sampling density compensation (CUDA specialization).
     * 
     * This class is intended for internal use only.
     * 
     * \tparam REAL Floating-point type.
     * \tparam D Number of dimensions.
     */
    template<class REAL, unsigned int D>
    class EXPORTSDC cuSDC_impl : public SDC_impl<cuNDArray, REAL, D>
    {

    public:

        cuSDC_impl(const vector_td<size_t, D>& matrix_size, REAL os_factor, size_t num_iterations, int device = -1);

    protected:
    
        virtual void preprocess(const cuNDArray<vector_td<REAL, D>>& traj) override;

        virtual void convolve(const cuNDArray<REAL>& in, cuNDArray<REAL>& out, SDC_conv_mode mode) override;

        virtual void update(const cuNDArray<REAL>& in, cuNDArray<REAL>& out) override;

    private:

        void convolve_NC2C(const cuNDArray<REAL>* in, cuNDArray<REAL>* out);
        void convolve_C2NC(const cuNDArray<REAL>* in, cuNDArray<REAL>* out);

        void check_consistency(const cuNDArray<REAL>* samples, const cuNDArray<REAL>* image, const cuNDArray<REAL>* dcw);

        void barebones();

        vector_td<size_t, D> grid_padding_;

        thrust::device_vector<vector_td<REAL, D>> traj_;

        int device_;
    };

    /**
     * \brief Estimate density compensation weights for an arbitrary trajectory.
     * 
     * \tparam REAL Floating-point type.
     * \tparam D Number of dimensions.
     * \param traj Flattened trajectory array of size [Ns, ...], where Ns is the number of samples in a single frame.
     * \param matrix_size Size of reconstructed image.
     * \param os_factor Oversampling factor
     * \param num_iterations Number of iterations.
     * \return boost::shared_ptr<cuNDArray<REAL>> Density compensation weights, same size as trajectory.
     */
    template<class REAL, unsigned int D>
    EXPORTSDC boost::shared_ptr<cuNDArray<REAL>> estimate_dcw(
        const cuNDArray<vector_td<REAL, D>>& traj,
        const vector_td<size_t, D>& matrix_size,
        REAL os_factor = 2.1,
        size_t num_iterations = 10);

    /**
     * \brief Estimate density compensation weights for an arbitrary trajectory.
     * 
     * \tparam REAL Floating-point type.
     * \tparam D Number of dimensions.
     * \param traj Flattened trajectory array of size [Ns, ...], where Ns is the number of samples in a single frame.
     * \param initial_dcw Initial estimate for density compensation weights.
     * \param matrix_size Size of reconstructed image.
     * \param os_factor Oversampling factor
     * \param num_iterations Number of iterations.
     * \return boost::shared_ptr<cuNDArray<REAL>> Density compensation weights, same size as trajectory.
     */
    template<class REAL, unsigned int D>
    EXPORTSDC boost::shared_ptr<cuNDArray<REAL>> estimate_dcw(
        const cuNDArray<vector_td<REAL, D>>& traj,
        const cuNDArray<REAL>& initial_dcw,
        const vector_td<size_t, D>& matrix_size,
        REAL os_factor = 2.1,
        size_t num_iterations = 10);

}   // namespace Gadgetron
