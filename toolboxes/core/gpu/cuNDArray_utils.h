/**
 * @file cuNDArray_utils.h
 */
#pragma once

#include "cuNDArray.h"
#include "vector_td.h"

namespace Gadgetron{

/**
 * @brief Cyclicly shifts the order of the array dimensions
 */
template<class T> cuNDArray<T>
shift_dim( const cuNDArray<T>& in, int shift );
/**
 * @brief Cyclicly shifts the order of the array dimensions
 */
template<class T> void
shift_dim( const cuNDArray<T>& in, cuNDArray<T>& out, int shift );

/**
 * @brief Permutes the array dimensions following the specified dimension order
 */
template<class T> cuNDArray<T>
permute(const cuNDArray<T>& in, const std::vector<size_t>& dim_order );

/**
 * @brief Permutes the array dimensions following the specified dimension order
 */
template<class T> void
permute( const cuNDArray<T>& in, cuNDArray<T>& out, const std::vector<size_t>& dim_order );

/**
 * @brief Creates a cropped version of the array
 * @param[in] crop_offset Offset of the corner of the crop size
 * @param[in] crop_size Size of the output array
 * @param[in] in Array to crop
 */
template<class T, unsigned int D> cuNDArray<T>
crop( const vector_td<size_t,D>& crop_offset, const vector_td<size_t,D>& crop_size, const cuNDArray<T>& in );

/**
 * @brief Creates a cropped version of the array
 * @param[in] crop_size Size of the cropped region
 * @param[in] in Array to crop
 * @param[out] out Array into which the cropped array is placed
 */
template<class T, unsigned int D> 
void crop( const vector_td<size_t,D>& crop_offset, const vector_td<size_t,D>& crop_size, const cuNDArray<T>& in, cuNDArray<T>& out );

/**
 * @brief Creates a padded version of the array
  * @param[in] size Size of the output array
 * @param[in] in Array to pad
 * @param[in] val Numerical value of the padding
 */
template<class T, unsigned int D> cuNDArray<T>
pad(const vector_td<size_t,D>& size, const cuNDArray<T>& in, T val = T(0) );


/**
 * @brief Creates a padded version of the array
 * @param[in] in Array to pad
 * @param[in] out Output array
 * @param[in] val Numerical value of the padding
 */
template<class T, unsigned int D> 
void pad( const cuNDArray<T>& in, cuNDArray<T>& out, T val = T(0) );

/**
 * @brief Fills the image with a given value outside a box
 * @param[in] matrix_size Box size
 * @param[in,out] image Array to fill
 * @param[in] val Fill value
 */
template<class T, unsigned int D> 
void fill_border( const vector_td<size_t,D>, cuNDArray<T>& image, T val = T(0) );

/**
 * @brief Fills the image with a given value outside a radius from the center
 * @param[in] radius Radius of the circle
 * @param[in,out] in_out Array to fill
 * @param[in] val Fill value
 */
template<class T, unsigned int D>
void fill_border( typename realType<T>::Type radius, cuNDArray<T>& in_out, T val= T(0) );

// Expand array to new dimension
/**
 * @brief Creates a new array, expanded into an additional dimension
 * @param[in] data Input data
 * @param[in] added_dim_size Size of the new dimension
 */
template<class T> cuNDArray<T>
expand(const cuNDArray<T>& data, size_t added_dim_size );

/**
 * @brief Creates an array of 2 times the size, created via linear interpolation
 * @param[in] in Array to upsample
 */
template<class T, unsigned int D> 
 cuNDArray<T> upsample( const cuNDArray<T>& in );

/**
 * @brief Creates an array of 2 times the size, created via linear interpolation
 * @param[in] in Array to upsample
 * @param[out] out Output array
 */
template<class T, unsigned int D> 
void upsample(const cuNDArray<T>& in, cuNDArray<T>& out );

/**
 * @brief Creates an array of half the size, created via linear interpolation
 * @param[in] in Array to downsample
 */
template<class T, unsigned int D> 
 cuNDArray<T> downsample( const cuNDArray<T>& in );

/**
 * @brief Creates an array of half the size, created via linear interpolation
 * @param[in] in Array to downsample
 * @param[out] out Output Array
 */
template<class T, unsigned int D> 
void downsample(const cuNDArray<T>& in, cuNDArray<T>& out );
}
