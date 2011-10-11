#pragma once
#include "gadgetron_export.h"

#include "cuNDArray.h"
#include "vector_td.h"
#include "vector_td_utilities.h"

#include <boost/smart_ptr.hpp>

//
// Utilities returning a shared_ptr to the resulting cuNDArray
// Component-wise operations.
//

EXPORTGPUCORE enum cuNDA_device { CUNDA_CURRENT_DEVICE, CUNDA_NDARRAY_DEVICE };

// Norm (float/double/complext arrays)
template<class REAL, class T> EXPORTGPUCORE 
boost::shared_ptr< cuNDArray<REAL> > 
cuNDA_norm( cuNDArray<T> *data, 
	    cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
	    cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Norm (reald arrays)
template<class REAL, unsigned int D> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<REAL> > 
cuNDA_norm( cuNDArray< typename reald<REAL,D>::Type > *data,
	    cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
	    cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Norm squared (float/double/complext arrays)
template<class REAL, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<REAL> > 
cuNDA_norm_squared( cuNDArray<T> *data,
		    cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
		    cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Norm squared (reald arrays)
template<class REAL, unsigned int D> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<REAL> > 
cuNDA_norm_squared( cuNDArray< typename reald<REAL,D>::Type > *data,
		    cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
		    cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Sum over dimension (scalar and vector_td arrays)
template<class T> EXPORTGPUCORE 
boost::shared_ptr< cuNDArray<T> >  
cuNDA_sum( cuNDArray<T> *data, unsigned int dim, 
	   cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
	   cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Sum over dimension (scalar and vector_td arrays)
template<class T> EXPORTGPUCORE 
boost::shared_ptr< cuNDArray<T> >  
cuNDA_expand( cuNDArray<T> *data, unsigned int added_dim_size, 
	      cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
	      cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Sum of Squares (float/double/complext array)
template<class S, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<S> > 
cuNDA_ss( cuNDArray<T> *data, unsigned int dim,
	  cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
	  cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Root Sum of Squares (float/double/complext array)
template<class S, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<S> > 
cuNDA_rss( cuNDArray<T> *data, unsigned int dim,
	   cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
	   cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Reciprocal Root Sum of Squares (float/double/complext array)
template<class S, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<S> > 
cuNDA_reciprocal_rss( cuNDArray<T> *data, unsigned int dim,
		      cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
		      cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Correlation matrix over the last dimension in the input array (float/double/complext array)
template<class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<T> > 
cuNDA_correlation( cuNDArray<T> *data,
		   cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
		   cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Real to complext
template<class REAL> EXPORTGPUCORE 
boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> > 
cuNDA_real_to_complext( cuNDArray<REAL> *data,
			cuNDA_device alloc_device = CUNDA_CURRENT_DEVICE, 
			cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

//
// Utilities overwriting the input
//

// Clear (real and complex types)
template<class T> EXPORTGPUCORE
bool cuNDA_clear( cuNDArray<T> *in_out, T val = get_zero<T>(), 
		  cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Reciprocal (real and complex types)
template<class T> EXPORTGPUCORE
bool cuNDA_reciprocal( cuNDArray<T> *in_out,
		       cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Square root (float/double for now)
template<class T> EXPORTGPUCORE
bool cuNDA_sqrt( cuNDArray<T> *in_out,
		 cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Reciprocal square root (float/double for now)
template<class T> EXPORTGPUCORE
bool cuNDA_reciprocal_sqrt( cuNDArray<T> *in_out,
			    cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Abs (floatd types supported, hence component-wise operation also for arrays of complex types!)
template<class T> EXPORTGPUCORE
bool cuNDA_abs( cuNDArray<T> *in_out,
		cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Normalize by RSS (float/double/complext arrays)
template<class REAL, class T> EXPORTGPUCORE
bool cuNDA_rss_normalize( cuNDArray<T> *in_out, unsigned int dim,
			  cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Scale (with constant)
template<class REAL> EXPORTGPUCORE
bool cuNDA_scale( REAL a, cuNDArray<typename complext<REAL>::Type> *x,
		  cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Scale (component-wise)
template<class T> EXPORTGPUCORE
bool cuNDA_scale( cuNDArray<T> *a, cuNDArray<T> *x,
		  cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

template<class T> EXPORTGPUCORE
bool cuNDA_scale_conj( cuNDArray<T> *a, cuNDArray<T> *x,
		  cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Scale (component-wise)
template<class REAL> EXPORTGPUCORE
bool cuNDA_scale( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x,
		  cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// 'axpy' - component-wise
template<class T> EXPORTGPUCORE
bool cuNDA_axpy( cuNDArray<T> *a, cuNDArray<T> *x, cuNDArray<T> *y,
		 cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// 'axpy' - component-wise
template<class REAL> EXPORTGPUCORE
bool cuNDA_axpy( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x, cuNDArray<typename complext<REAL>::Type> *y,
		 cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

//
// cublas wrappers (real and complex types)
//

template<class T> EXPORTGPUCORE 
T cuNDA_dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2,
	     cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

template<class REAL, class T> EXPORTGPUCORE
REAL cuNDA_asum( cuNDArray<T>* arr,
		 cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

template<class T> EXPORTGPUCORE
bool cuNDA_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y,
		 cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

template<class T> EXPORTGPUCORE
bool cuNDA_scal( T a, cuNDArray<T>* x,
		 cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

template<class T> EXPORTGPUCORE
bool cuNDA_add( T a, cuNDArray<T>* x,
		 cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );


// Normalize (float/double arrays only)
template<class REAL> EXPORTGPUCORE
REAL cuNDA_normalize( cuNDArray<REAL> *in_out, REAL new_max,
		      cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

//
// Some image utilities (with an interface fitting the NNFT)
//

// Crop (scalar and vector_td arrays)
template<class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_crop( typename uintd<D>::Type crop_offset, 
		 cuNDArray<T> *in, cuNDArray<T> *out,
		 cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Expand with zero filling (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out,
				  cuNDA_device compute_device = CUNDA_CURRENT_DEVICE );

// Zero fill border (rectangular) - (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_zero_fill_border( typename uintd<D>::Type matrix_size, cuNDArray<T> *image,
			     cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Border fill (circular) - (real and complex types)
template<class REAL, class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_zero_fill_border( typename reald<REAL,D>::Type radius, cuNDArray<T> *image,
			     cuNDA_device compute_device = CUNDA_NDARRAY_DEVICE );

// Shrinkage operators
template<class REAL, class T> EXPORTGPUCORE
bool cuNDA_shrink1( REAL gamma, cuNDArray<T> *in, cuNDArray<T> *out );

template<class REAL, class T> EXPORTGPUCORE
bool cuNDA_shrinkd( REAL gamma, cuNDArray<REAL> *s_k, cuNDArray<T> *in, cuNDArray<T> *out );
