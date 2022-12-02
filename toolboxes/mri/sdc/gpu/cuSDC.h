#pragma once

#include "SDC.h"
#include "cuGriddingConvolution.h"
namespace Gadgetron {

/**
 * \brief Estimate density compensation weights for an arbitrary trajectory.
 *
 * \tparam ARRAY Array type.
 * \tparam REAL Floating-point type.
 * \tparam D Number of dimensions.
 * \param traj Flattened trajectory array of size [Ns, ...], where Ns is the number of samples in a single frame.
 * \param matrix_size Size of reconstructed image.
 * \param os_factor Oversampling factor
 * \param num_iterations Number of iterations.
 * \param convtype is ConvolutionType::STANDARD ConvolutionType::ATOMIC
 * \return std::shared_ptr<ARRAY<REAL>> Density compensation weights, same size as trajectory.
 */
template <class REAL, unsigned int D>
std::shared_ptr<cuNDArray<REAL>>
estimate_dcw(const cuNDArray<vector_td<REAL, D>>& traj, const vector_td<size_t, D>& matrix_size, REAL os_factor = 2.1,
             unsigned int num_iterations = 10, ConvolutionType convtype = ConvolutionType::STANDARD);

/**
 * \brief Estimate density compensation weights for an arbitrary trajectory.
 *
 * \tparam ARRAY Array type.
 * \tparam REAL Floating-point type.
 * \tparam D Number of dimensions.
 * \param traj Flattened trajectory array of size [Ns, ...], where Ns is the number of samples in a single frame.
 * \param initial_dcw Initial estimate for density compensation weights.
 * \param matrix_size Size of reconstructed image.
 * \param os_factor Oversampling factor
 * \param num_iterations Number of iterations.
 * \param convtype is ConvolutionType::STANDARD ConvolutionType::ATOMIC
 * \return std::shared_ptr<ARRAY<REAL>> Density compensation weights, same size as trajectory.
 */
template <class REAL, unsigned int D>
std::shared_ptr<cuNDArray<REAL>>
estimate_dcw(const cuNDArray<vector_td<REAL, D>>& traj, const cuNDArray<REAL>& initial_dcw,
             const vector_td<size_t, D>& matrix_size, REAL os_factor = 2.1, unsigned int num_iterations = 10,
             ConvolutionType convtype = ConvolutionType::STANDARD);
} // namespace Gadgetron