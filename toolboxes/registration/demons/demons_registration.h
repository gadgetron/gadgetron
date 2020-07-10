//
// Created by dchansen on 1/29/20.
//

#pragma once
#include "hoNDArray.h"

namespace Gadgetron {
namespace Registration {

/**
 * This function takes an image and deforms it by the vector field using linear
 * interpolation
 * @tparam T datatype of the image. Must be floating point like, such as double,
 * vector_td<float,3> or std::complex<float>
 * @tparam R Datatype of the deformation field. Must be floating point (float or
 * double)
 * @tparam D Dimensionality of the images
 * @param image Image to deform
 * @param deformation_field Deformation field by which each voxel should be
 * offset.
 * @return The deformed image.
 */
template <class T, unsigned int D, class R = realType_t<T>>
hoNDArray<T> deform_image(const hoNDArray<T> &image,
                          const hoNDArray<vector_td<R, D>> &deformation_field);

template <class T, unsigned int D, class R = realType_t<T>>
hoNDArray<T>
deform_image_bspline(const hoNDArray<T> &image,
                     const hoNDArray<vector_td<R, D>> &deformation_field);
/**
 *
 * @tparam T datatype of the image. Float or double
 * @tparam D Image dimension.
 * @param fixed
 * @param moving
 * @param iterations Number of iterations
 * @param sigma Sigma for smoothing the motion field
 * @return The vector field deforming fixed to moving.
 */
template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
diffeomorphic_demons(const hoNDArray<T> &fixed, const hoNDArray<T> &moving,
                     unsigned int iterations = 20, float sigma = 2.0,
                     float step_size = 2.0, float noise_sigma = 0.0f);

/**
 *
 * @tparam T datatype of the image. Float or double
 * @tparam D Image dimension.
 * @param vector_field Deformation field used as initial guess.
 * @param fixed
 * @param moving
 * @param iterations Number of iterations
 * @param sigma Sigma for smoothing the motion field
 * @return The vector field deforming fixed to moving.
 */
template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
diffeomorphic_demons(const hoNDArray<T> &fixed, const hoNDArray<T> &moving,
                     hoNDArray<vector_td<T, D>> vector_field,
                     unsigned int iterations = 20, float sigma = 2.0,
                     float step_size = 2.0, float noise_sigma = 0.0f);

/**
 *
 * @tparam T datatype of the image. Float or double
 * @tparam D Image dimension.
 * @param fixed
 * @param moving
 * @param iterations Number of iterations
 * @param sigma Sigma for smoothing the motion field
 * @param step_size Step size for every iteration
 * @param gradient_eps
 * @return The vector field deforming fixed to moving.
 */
template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
ngf_diffeomorphic_demons(const hoNDArray<T> &fixed, const hoNDArray<T> &moving,
                         unsigned int iterations = 20, float sigma = 2.0,
                         float step_size = 2, float gradient_eps = 1e-6f);


/**
 *
 * @tparam T datatype of the image. Float or double
 * @tparam D Image dimension.
 * @param vector_field Deformation field used as initial guess.
 * @param fixed
 * @param moving
 * @param iterations Number of iterations
 * @param sigma Sigma for smoothing the motion field
 * @param step_size Step size for every iteration
 * @param gradient_eps
 * @return The vector field deforming fixed to moving.
 */
template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
ngf_diffeomorphic_demons(const hoNDArray<T> &fixed, const hoNDArray<T> &moving,hoNDArray<vector_td<T,D>> vector_field,
                         unsigned int iterations = 20, float sigma = 2.0,
                         float step_size = 2, float gradient_eps = 1e-6f);

/**
 *
 * @tparam T datatype of the image. Float or double
 * @tparam D Image dimension.
 * @param vector_field Deformation field used as initial guess.
 * @param fixed
 * @param moving
 * @param levels Number of pyramid levels
 * @param iterations Number of iterations
 * @param sigma Sigma for smoothing the motion field
 * @return The vector field deforming fixed to moving.
 */
template <class T, unsigned int D>
hoNDArray<vector_td<T, D>> multi_scale_diffeomorphic_demons(
    const hoNDArray<T> &fixed, const hoNDArray<T> &moving,
    unsigned int levels = 3, unsigned int iterations = 20, float sigma = 2.0,
    float step_size = 2.0, float noise_sigma = 0.0);


/**
 *
 * @tparam T datatype of the image. Float or double
 * @tparam D Image dimension.
 * @param vector_field Deformation field used as initial guess.
 * @param fixed
 * @param moving
 * @param levels Number of pyramid levels
 * @param iterations Number of iterations
 * @param sigma Sigma for smoothing the motion field
 * @param step_size Step sized used in the demons algorithm
 * @param gradient_eps Constant to add when calculating the NGF denominator
 * @return The vector field deforming fixed to moving.
 */
template <class T, unsigned int D>
hoNDArray<vector_td<T, D>> multi_scale_ngf_diffeomorphic_demons(
    const hoNDArray<T> &fixed, const hoNDArray<T> &moving,
    unsigned int levels = 3, unsigned int iterations = 20, float sigma = 2.0,
    float step_size = 2.0, float gradient_eps = 1e-6);
template <class T>
hoNDArray<T> gaussian_filter(const hoNDArray<T> &image, float sigma);

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
compose_fields(const hoNDArray<vector_td<T, D>> &update_field,
               const hoNDArray<vector_td<T, D>> &vfield);

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
vector_field_exponential(const hoNDArray<vector_td<T, D>> &vector_field);

hoNDArray<vector_td<float, 2>> demons_step_ext(const hoNDArray<float> &fixed,
                                               const hoNDArray<float> &moving,
                                               float alpha, float beta,
                                               float noise_sigma);
} // namespace Registration
} // namespace Gadgetron
