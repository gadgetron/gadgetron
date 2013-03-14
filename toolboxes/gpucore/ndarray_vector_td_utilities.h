#pragma once
#include "gpucore_export.h"

#include "cuNDArray.h"
#include "vector_td.h"
#include "vector_td_utilities.h"
#include "complext.h"
#include <boost/smart_ptr.hpp>

//
// Utilities returning a shared_ptr to the resulting cuNDArray
// Component-wise operations.
//

enum cuNDA_device {
	CUNDA_CURRENT_DEVICE, CUNDA_NDARRAY_DEVICE
};


namespace Gadgetron{
// Abs "complex style" (float/double/complext arrays)
/**
 * @brief Calculates the elementwise absolute value of the array
 * @param[in] data Input data
 * @param[in] alloc_device Device on which to allocate the new array
 * @param[in] compute_device Device on which to do the computation
 * @return A new array containing the elementwise absolute value of data
 */
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<typename realType<T>::type> >
abs(cuNDArray<T> *data, cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Sum over dimension (scalar and vector_td arrays)
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<T> >
sum(cuNDArray<T> *data, unsigned int dim, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Expand (copy) array to new dimension (scalar and vector_td arrays)
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<T> >
expand(cuNDArray<T> *data, unsigned int added_dim_size,
		cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);


// Correlation matrix over the last dimension in the input array (float/double/complext array)
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<T> >
correlation(cuNDArray<T> *data,
		cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Real to complext
template<class REAL> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<complext<REAL> > >
real_to_complext(cuNDArray<REAL> *data, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// complext to real (by discarding the imaginary component)
template<class REAL> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<REAL> >
complext_to_real(cuNDArray<complext<REAL> > *data, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Downsample array to half size (Real arrays only)
template<class REAL, unsigned int D> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<REAL> >
downsample(cuNDArray<REAL> *data, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Nearest neighbor upsampling of array to double size (Real arrays only)
template<class REAL, unsigned int D> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<REAL> >
upsample_nn(cuNDArray<REAL> *data, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Linear interpolation upsampling of array to double size (Real arrays only)
template<class REAL, unsigned int D> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<REAL> >
upsample_lin(cuNDArray<REAL> *data, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

/**
 * @brief Clamps all values in the array to the minimum and maximum values specified.
 * @param[in,out] in_out Array which to clamp
 * @param[in] min minimum value
 * @param[in] max maximum value
 */
template<class  T> EXPORTGPUCORE
void clamp(cuNDArray<T> *in_out, T min, T max);


/**
 * @brief Clamps all values in the array to the minimum value specified.
 * @param[in,out] in_out Array which to clamp
 * @param[in] min minimum value
 */
template<class  T> EXPORTGPUCORE
void clamp_min(cuNDArray<T> *in_out, T min);

/**
 * @brief Clamps all values in the array to the maximum value specified.
 * @param[in,out] in_out Array which to clamp
 * @param[in] max minimum value
 */
template<class  T> EXPORTGPUCORE
void clamp_max(cuNDArray<T> *in_out, T max);

// Normalize by RSS (float/double/complext arrays)
template<class T> EXPORTGPUCORE
void rss_normalize(cuNDArray<T> *in_out, unsigned int dim,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

//Elementwise sgn function.
/**
 * @brief Calculates the signum function elementwise
 * @param[in,out] data
 */
template <class T> EXPORTGPUCORE void inplace_sgn(cuNDArray<T>* data);


// Sum of Squares (float/double/complext array)
/**
 * @brief Calculates the squared norm along the specified dimension
 * @param[in] data Input data
 * @param[in] dim Dimension along which to do the squared norm
 * @param[in] alloc_device Device on which to allocate the return array
 * @param[in] compute_device Device on which to do the computations
 * @return boost shared pointer to array containing the elementwise squared norm of data along the specifie dimension
 */
template<class T>  EXPORTGPUCORE boost::shared_ptr<cuNDArray<typename realType<T>::type > >
squaredNorm(cuNDArray<T> *data, unsigned int dim, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);
//
// Some image utilities (with an interface fitting the NNFT)
//
/**
 * Calculates the elementwise maximum of two arrays
 * @param[in] in1 First input array
 * @param[in] in2 Second input Array
 * @param[in] alloc_device Device on which to allocate the new array
 * @param[in] compute_device Device on which to do the computation
 * @return shared pointer to array containing the elementwise maximum of two arrays
 */
template<class T>
boost::shared_ptr< cuNDArray<T> >
maximum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device );
/**
 * Calculates the elementwise minimum of two arrays
 * @param[in] in1 First input array
 * @param[in] in2 Second input Array
 * @param[in] alloc_device Device on which to allocate the new array
 * @param[in] compute_device Device on which to do the computation
 * @return shared pointer to array containing the elementwise minimum of two arrays
 */
template<class T>
boost::shared_ptr< cuNDArray<T> >
minimum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device );
// Crop (scalar and vector_td arrays)
template<class T, unsigned int D> EXPORTGPUCORE
void crop(typename uintd<D>::Type crop_offset, cuNDArray<T> *in,
		cuNDArray<T> *out, cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Expand with zero filling (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
void expand_with_zero_fill(cuNDArray<T> *in, cuNDArray<T> *out,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Zero fill border (rectangular) - (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
void zero_fill_border(typename uintd<D>::Type matrix_size, cuNDArray<T> *image,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Border fill (circular) - (real and complex types)
template<class REAL, class T, unsigned int D> EXPORTGPUCORE
void zero_fill_border(REAL radius, cuNDArray<T> *image,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Shrinkage operators
template<class REAL, class T> EXPORTGPUCORE
void shrink1(REAL gamma, cuNDArray<T> *in, cuNDArray<T> *out);

template<class REAL, class T> EXPORTGPUCORE
void shrinkd(REAL gamma, cuNDArray<REAL> *s_k, cuNDArray<T> *in,
		cuNDArray<T> *out);

// Mirror around the origin -- !! leaving the origin unchanged !!
template<class T, unsigned int D> EXPORTGPUCORE
void origin_mirror(cuNDArray<T> *in, cuNDArray<T> *out, bool zero_fill = true,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);


template<class T> EXPORTGPUCORE
T normalize( cuNDArray<T> *data, T new_max, cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

}
