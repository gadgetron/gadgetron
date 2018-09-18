/**
 * @file cuNDArray_utils.h
 */
#pragma once

#include "cuNDArray.h"
#include "vector_td.h"
#include "gpucore_export.h"

namespace Gadgetron{

/**
 * @brief Cyclicly shifts the order of the array dimensions
 */
template<class T> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
shift_dim( cuNDArray<T> *in, int shift );
/**
 * @brief Cyclicly shifts the order of the array dimensions
 */
template<class T> EXPORTGPUCORE void
shift_dim( cuNDArray<T> *in, cuNDArray<T> *out, int shift );

/**
 * @brief Permutes the array dimensions following the specified dimension order
 */
template<class T> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
permute( cuNDArray<T> *in, std::vector<size_t> *dim_order, int shift_mode = 0 );

/**
 * @brief Permutes the array dimensions following the specified dimension order
 */
template<class T> EXPORTGPUCORE void
permute( cuNDArray<T> *in, cuNDArray<T> *out, std::vector<size_t> *dim_order, int shift_mode = 0 );

/**
 * @brief Creates a cropped version of the array
 * @param[in] crop_offset Offset of the corner of the crop size
 * @param[in] crop_size Size of the output array
 * @param[in] in Array to crop
 */
template<class T, unsigned int D> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
crop( typename uint64d<D>::Type crop_offset, typename uint64d<D>::Type crop_size, cuNDArray<T> *in );

/**
 * @brief Creates a cropped version of the array
 * @param[in] crop_offset Offset of the corner of the crop size
 * @param[in] in Array to crop
 * @param[out] out Array into which the cropped array is placed
 */
template<class T, unsigned int D> EXPORTGPUCORE
void crop( typename uint64d<D>::Type crop_offset, cuNDArray<T> *in, cuNDArray<T> *out );

/**
 * @brief Creates a padded version of the array
  * @param[in] size Size of the output array
 * @param[in] in Array to pad
 * @param[in] val Numerical value of the padding
 */
template<class T, unsigned int D> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
pad( typename uint64d<D>::Type size, cuNDArray<T> *in, T val = T(0) );


/**
 * @brief Creates a padded version of the array
 * @param[in] in Array to pad
 * @param[in] out Output array
 * @param[in] val Numerical value of the padding
 */
template<class T, unsigned int D> EXPORTGPUCORE
void pad( cuNDArray<T> *in, cuNDArray<T> *out, T val = T(0) );

/**
 * @brief Fills the image with a given value outside a box
 * @param[in] matrix_size Box size
 * @param[in,out] image Array to fill
 * @param[in] val Fill value
 */
template<class T, unsigned int D> EXPORTGPUCORE
void fill_border( typename uint64d<D>::Type matrix_size, cuNDArray<T> *image, T val = T(0) );

/**
 * @brief Fills the image with a given value outside a radius from the center
 * @param[in] radius Radius of the circle
 * @param[in,out] in_out Array to fill
 * @param[in] val Fill value
 */
template<class T, unsigned int D>
void fill_border( typename realType<T>::Type radius, cuNDArray<T> *in_out, T val= T(0) );

// Expand array to new dimension
/**
 * @brief Creates a new array, expanded into an additional dimension
 * @param[in] data Input data
 * @param[in] added_dim_size Size of the new dimension
 */
template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> >
expand(cuNDArray<T> *data, size_t added_dim_size );

/**
 * @brief Creates an array of 2 times the size, created via linear interpolation
 * @param[in] in Array to upsample
 */
template<class T, unsigned int D> EXPORTGPUCORE
 cuNDArray<T> upsample( cuNDArray<T>* in );

/**
 * @brief Creates an array of 2 times the size, created via linear interpolation
 * @param[in] in Array to upsample
 * @param[out] out Output array
 */
template<class T, unsigned int D> EXPORTGPUCORE
void upsample( cuNDArray<T> *in, cuNDArray<T> *out );

/**
 * @brief Creates an array of half the size, created via linear interpolation
 * @param[in] in Array to downsample
 */
template<class T, unsigned int D> EXPORTGPUCORE
 cuNDArray<T> downsample( cuNDArray<T>* in );

/**
 * @brief Creates an array of half the size, created via linear interpolation
 * @param[in] in Array to downsample
 * @param[out] out Output Array
 */
template<class T, unsigned int D> EXPORTGPUCORE
void downsample( cuNDArray<T> *in, cuNDArray<T> *out );
}
