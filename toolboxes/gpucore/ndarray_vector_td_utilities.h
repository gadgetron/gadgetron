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

// Abs "complex style" (float/double/complext arrays)
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
/*
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray< typename realType<T>::type > >
norm(cuNDArray<T> *data, unsigned int dim, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

*/
// Sum of Squares (float/double/complext array)
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<typename realType<T>::type > >
squaredNorm(cuNDArray<T> *data, unsigned int dim, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);
/*
// Root Sum of Squares (float/double/complext array)
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray<typename realType<T>::type> >
rss(cuNDArray<T> *data, unsigned int dim, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);

// Reciprocal Root Sum of Squares (float/double/complext array)
 * */
/*
template<class T> EXPORTGPUCORE
boost::shared_ptr<cuNDArray< typename realType<T>::type> >
reciprocal_rss(cuNDArray<T> *data, unsigned int dim, cuNDA_device alloc_device =
		CUNDA_CURRENT_DEVICE,
		cuNDA_device compute_device = CUNDA_CURRENT_DEVICE);
*/
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

//
// Utilities overwriting the input
//

// Reciprocal square root (float/double for now)
/*template<class T> EXPORTGPUCORE
void reciprocal_sqrt(cuNDArray<T> *in_out, cuNDA_device compute_device =
		CUNDA_NDARRAY_DEVICE);*/


// Threshold
template<class T> EXPORTGPUCORE
void threshold_min(T min, cuNDArray<T> *in_out, T to_value = T(0),
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Threshold
template<class T> EXPORTGPUCORE
void threshold_max(T max, cuNDArray<T> *in_out, T to_value = T(0),
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Threshold
template<class T> EXPORTGPUCORE
void threshold_min(cuNDArray<T> * min, cuNDArray<T> *in_out, T to_value = T(0),
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Threshold
template<class T> EXPORTGPUCORE
void threshold_max(cuNDArray<T> * max, cuNDArray<T> *in_out, T to_value = T(0),
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);
template<class T> EXPORTGPUCORE
void threshold_amin(cuNDArray<T> * min, cuNDArray<T> *in_out, T to_value = T(0),
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

// Normalize by RSS (float/double/complext arrays)
template<class T> EXPORTGPUCORE
void rss_normalize(cuNDArray<T> *in_out, unsigned int dim,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

/*
// Scale (with constant)
template<class REAL> EXPORTGPUCORE
void scal(REAL a, cuNDArray<complext<REAL> > *x, cuNDA_device compute_device =
		CUNDA_NDARRAY_DEVICE);

// Scale (component-wise)
template<class T> EXPORTGPUCORE
void scale(cuNDArray<T> *a, cuNDArray<T> *x, cuNDA_device compute_device =
		CUNDA_NDARRAY_DEVICE);

template<class T> EXPORTGPUCORE
void scale_conj(cuNDArray<T> *a, cuNDArray<T> *x, cuNDA_device compute_device =
		CUNDA_NDARRAY_DEVICE);

// Scale (component-wise)
template<class REAL> EXPORTGPUCORE
void scale(cuNDArray<REAL> *a, cuNDArray<complext<REAL> > *x,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE);

*/

//
// Some image utilities (with an interface fitting the NNFT)
//

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
