#pragma once
#include "gadgetron_export.h"

#include "cuNDArray.h"
#include "vector_td.h"

#include <vector>

#include <cublas_v2.h>
#include <boost/smart_ptr.hpp>

//
// Utilities returning a shared_ptr to the resulting cuNDArray
// Component-wise operations.
//

// Sum over dimension dim (scalar and vector_td arrays)
template<class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<T> > cuNDA_sum( cuNDArray<T> *data, unsigned int dim );

// Norm (float/double/complext arrays)
template<class REAL, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<T> *data );

// Norm (reald arrays)
template<class REAL, unsigned int D> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray< typename reald<REAL,D>::Type > *data );

// Norm squared (float/double/complext arrays)
template<class REAL, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<T> *data );

// Norm squared (reald arrays)
template<class REAL, unsigned int D> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray< typename reald<REAL,D>::Type > *data );

// Sum of Squares (float/double/complext array)
template<class S, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<S> > cuNDA_ss( cuNDArray<T> *data, unsigned int dim );

// Root Sum of Squares (float/double/complext array)
template<class S, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<S> > cuNDA_rss( cuNDArray<T> *data, unsigned int dim );

// Reciprocal RSS (float/double/complext array)
template<class S, class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<S> > cuNDA_reciprocal_rss( cuNDArray<T> *data, unsigned int dim );

// Correlation matrix (float/double/complext array)
template<class T> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<T> > cuNDA_correlation( cuNDArray<T> *data );

// Real to complext
template<class REAL> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> > cuNDA_real_to_complext( cuNDArray<REAL> *data );

//
// Utilities overwriting the input
//

// Clear (real and complex types)
template<class T> EXPORTGPUCORE
void cuNDA_clear( cuNDArray<T> *in_out );

// Reciprocal (real and complex types)
template<class T> EXPORTGPUCORE
void cuNDA_reciprocal( cuNDArray<T> *in_out );

// Square root (float/double for now)
template<class T> EXPORTGPUCORE
void cuNDA_sqrt( cuNDArray<T> *in_out );

// Reciprocal square root (float/double for now)
template<class T> EXPORTGPUCORE
void cuNDA_reciprocal_sqrt( cuNDArray<T> *in_out );

// Normalize (float/double arrays)
template<class REAL> EXPORTGPUCORE
void cuNDA_normalize( cuNDArray<REAL> *in_out, REAL new_max, cublasHandle_t handle );

// Abs (floatd types supported, hence component-wise operation also for arrays of complex types!)
template<class T> EXPORTGPUCORE
void cuNDA_abs( cuNDArray<T> *in_out );

// Normalize by RSS (float/double/complext arrays)
template<class REAL, class T> EXPORTGPUCORE
bool cuNDA_rss_normalize( cuNDArray<T> *in_out, unsigned int dim );

// Scale (with constant)
template<class REAL> EXPORTGPUCORE
void cuNDA_scale( REAL a, cuNDArray<typename complext<REAL>::Type> *x );

// Scale (component-wise)
template<class T> EXPORTGPUCORE
bool cuNDA_scale( cuNDArray<T> *a, cuNDArray<T> *x );

// Scale (component-wise)
template<class REAL> EXPORTGPUCORE
bool cuNDA_scale( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x );

// 'axpy' - component-wise
template<class T> EXPORTGPUCORE
bool cuNDA_axpy( cuNDArray<T> *a, cuNDArray<T> *x, cuNDArray<T> *y );

// 'axpy' - component-wise
template<class REAL> EXPORTGPUCORE
bool cuNDA_axpy( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x, cuNDArray<typename complext<REAL>::Type> *y );

// Crop
template<class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_crop( typename uintd<D>::Type offset, cuNDArray<T> *in, cuNDArray<T> *out );

// Expand with zero filling (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out );

// Zero fill border (rectangular) - (real and complex types)
template<class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_zero_fill_border( typename uintd<D>::Type matrix_size, cuNDArray<T> *image );

// Border fill (circular) - (real and complex types)
template<class REAL, class T, unsigned int D> EXPORTGPUCORE
bool cuNDA_zero_fill_border( typename reald<REAL,D>::Type radius, cuNDArray<T> *image );

//
// cublas wrappers (real and complex types)
//

template<class T> EXPORTGPUCORE 
T cuNDA_dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, cublasHandle_t handle );

template<class REAL, class T> EXPORTGPUCORE
REAL cuNDA_asum( cuNDArray<T>* arr, cublasHandle_t handle );

template<class T> EXPORTGPUCORE
bool cuNDA_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y, cublasHandle_t handle );

template<class T> EXPORTGPUCORE
bool cuNDA_scal( T a, cuNDArray<T>* x, cublasHandle_t handle );
