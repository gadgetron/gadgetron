#pragma once

#include "cuNDArray.h"
#include "vector_td.h"

#include <memory>
#include <vector>

#include <cublas_v2.h>

//
// Utilities returning an auto_ptr to the resulting cuNDArray
// Component-wise operations.
//

// Sum over dimension dim (scalar and vector_td arrays)
template<class T>
std::auto_ptr< cuNDArray<T> > cuNDA_sum( cuNDArray<T> *data, unsigned int dim );

// Norm (float/double/complext arrays)
template<class REAL, class T>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<T> *data );

// Norm (reald arrays)
template<class REAL, unsigned int D>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray< typename reald<REAL,D>::Type > *data );

// Norm squared (float/double/complext arrays)
template<class REAL, class T>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<T> *data );

// Norm squared (reald arrays)
template<class REAL, unsigned int D>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray< typename reald<REAL,D>::Type > *data );

// RSS (float/double/complext array)
template<class S, class T>
std::auto_ptr< cuNDArray<S> > cuNDA_rss( cuNDArray<T> *data, unsigned int dim );

// Reciprocal RSS (float/double/complext array)
template<class S, class T>
std::auto_ptr< cuNDArray<S> > cuNDA_reciprocal_rss( cuNDArray<T> *data, unsigned int dim );

// Correlation matrix (float/double/complext array)
template<class T>
std::auto_ptr< cuNDArray<T> > cuNDA_correlation( cuNDArray<T> *data );

//
// Utilities overwriting the input
//

// Clear (real and complex types)
template<class T> 
void cuNDA_clear( cuNDArray<T> *in_out );

// Reciprocal (real and complex types)
template<class T> 
void cuNDA_reciprocal( cuNDArray<T> *in_out );

// Normalize (float/double arrays)
template<class REAL>
void cuNDA_normalize( cuNDArray<REAL> *in_out, REAL new_max, cublasHandle_t handle );

// Abs (floatd types supported, hence component-wise operation also for arrays of complex types!)
template<class T> 
void cuNDA_abs( cuNDArray<T> *in_out );

// Normalize by RSS (float/double/complext arrays)
template<class REAL, class T> 
bool cuNDA_rss_normalize( cuNDArray<T> *in_out, unsigned int dim );

// Scale (with constant)
template<class REAL> 
void cuNDA_scale( REAL a, cuNDArray<typename complext<REAL>::Type> *x );

// Scale (component-wise)
template<class T>
bool cuNDA_scale( cuNDArray<T> *a, cuNDArray<T> *x );

// Scale (component-wise)
template<class REAL>
bool cuNDA_scale( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x );

// 'axpby' - overwrites y
template<class T>
bool cuNDA_axpby( cuNDArray<T> *a, cuNDArray<T> *x, cuNDArray<T> *b, cuNDArray<T> *y );

// 'axpby' - overwrites y
template<class REAL>
bool cuNDA_axpby( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x, cuNDArray<REAL> *b, cuNDArray<typename complext<REAL>::Type> *y );

// Crop
template<class T, unsigned int D>
bool cuNDA_crop( typename uintd<D>::Type offset, cuNDArray<T> *in, cuNDArray<T> *out );

// Expand with zero filling (real and complex types)
template<class T, unsigned int D>
bool cuNDA_expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out );

// Zero fill border (rectangular) - (real and complex types)
template<class T, unsigned int D> 
bool cuNDA_zero_fill_border( typename uintd<D>::Type matrix_size, cuNDArray<T> *image );

// Border fill (circular) - (real and complex types)
template<class REAL, class T, unsigned int D>
bool cuNDA_zero_fill_border( typename reald<REAL,D>::Type radius, cuNDArray<T> *image );

//
// cublas wrappers (real and complex types)
//

template<class T> T
cuNDA_dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, cublasHandle_t handle );

template<class T> bool
cuNDA_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y, cublasHandle_t handle );

template<class T> bool
cuNDA_scal( T a, cuNDArray<T>* x, cublasHandle_t handle );
