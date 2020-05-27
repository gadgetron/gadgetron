/**
 * \file hoSDC.h
 * \brief Sampling density compensation (host specialization).
 */

#pragma once

#include "SDC.h"
#include "mri_sdc_export.h"

#include "hoSDC_conv.h"

namespace Gadgetron
{
    /**
     * \brief Implementation of sampling density compensation (host specialization).
     * 
     * This class is intended for internal use only.
     * 
     * \tparam REAL Floating-point type.
     * \tparam D Number of dimensions.
     */
    template<class REAL, unsigned int D>
    class EXPORTSDC hoSDC_impl : public SDC_impl<hoNDArray, REAL, D>
    {
    public:

        hoSDC_impl(const vector_td<size_t, D>& matrix_size, REAL os_factor, size_t num_iterations);

    protected:

        virtual void preprocess(const hoNDArray<vector_td<REAL, D>>& traj) override;

        virtual void convolve(const hoNDArray<REAL>& in, hoNDArray<REAL>& out, SDC_conv_mode mode) override;

        virtual void update(const hoNDArray<REAL>& in, hoNDArray<REAL>& out) override;

    private:

        void convolve_NC2C(const hoNDArray<REAL>& in, hoNDArray<REAL>& out);
        void convolve_C2NC(const hoNDArray<REAL>& in, hoNDArray<REAL>& out);

        std::vector<SDC_internal::hoConvMatrix<REAL>> conv_matrix_;
        std::vector<SDC_internal::hoConvMatrix<REAL>> conv_matrix_T_;
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
     * \return boost::shared_ptr<hoNDArray<REAL>> Density compensation weights, same size as trajectory.
     */
    template<class REAL, unsigned int D>
    EXPORTSDC boost::shared_ptr<hoNDArray<REAL>> estimate_dcw(
        const hoNDArray<vector_td<REAL, D>>& traj,
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
     * \return boost::shared_ptr<hoNDArray<REAL>> Density compensation weights, same size as trajectory.
     */
    template<class REAL, unsigned int D>
    EXPORTSDC boost::shared_ptr<hoNDArray<REAL>> estimate_dcw(
        const hoNDArray<vector_td<REAL, D>>& traj,
        const hoNDArray<REAL>& initial_dcw,
        const vector_td<size_t, D>& matrix_size,
        REAL os_factor = 2.1,
        size_t num_iterations = 10);

}   // namespace Gadgetron
