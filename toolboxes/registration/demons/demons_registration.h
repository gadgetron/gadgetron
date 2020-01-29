//
// Created by dchansen on 1/29/20.
//

#pragma once
#include <hoNDArray.h>

namespace Gadgetron { namespace Registration {

/**
 * This function takes an image and deforms it by the vector field using linear interpolation
 * @tparam T datatype of the image. Must be floating point like, such as double, vector_td<float,3> or std::complex<float>
 * @tparam R Datatype of the deformation field. Must be floating point (float or double)
 * @tparam D Dimensionality of the images
 * @param image Image to deform
 * @param deformation_field Deformation field by which each voxel should be offset.
 * @return The deformed image.
 */
    template <class T, class R, unsigned int D>
    hoNDArray<T> deform_image(const hoNDArray<T>& image, const hoNDArray<vector_td<R, D>>& deformation_field);


}}
