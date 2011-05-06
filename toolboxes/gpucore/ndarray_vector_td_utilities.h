#pragma once

#include "cuNDArray.h"
#include "vector_td.h"

#include <memory>
#include <vector>

#include <cublas_v2.h>

//	
// cuNDArray of scalar/vector_td utilities
//

//
// Utilities returning an auto_ptr to the resulting cuNDArray
//

// Sum (scalar and vector_td types)
template<class T>
std::auto_ptr< cuNDArray<T> > cuNDA_sum( cuNDArray<T> *data, unsigned int dim );

// Norm (real types)
template<class REAL>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<REAL> *data );

// Norm (reald types)
template<class REAL, unsigned int D>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<typename reald<REAL,D>::Type> *data );

// Norm squared (real types)
template<class REAL>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<REAL> *data );

// Norm squared (reald types)
template<class REAL, unsigned int D>
std::auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<typename reald<REAL,D>::Type> *data );

// RSS (real and complext types)
template<class REAL, class T>
std::auto_ptr< cuNDArray<REAL> > cuNDA_rss( cuNDArray<T> *data, unsigned int dim );

// Correlation matrix (complext types)
template<class REAL>
std::auto_ptr< cuNDArray<typename complext<REAL>::Type> > cuNDA_correlation( cuNDArray<typename complext<REAL>::Type> *data );

//
// Utilities overwriting the input
//

// Clear (real and complex types)
template<class T> 
void cuNDA_clear( cuNDArray<T> *in_out );

// Reciprocal (real and complex types)
template<class T> 
void cuNDA_reciprocal( cuNDArray<T> *in_out );

// Normalize (real types)
// TODO: complex version
template<class REAL>
void cuNDA_normalize( cuNDArray<REAL> *in_out, REAL new_max, cublasHandle_t handle );

// Abs (floatd types supported, hence component-wise operation also for arrays of complex types!)
template<class T> 
void cuNDA_abs( cuNDArray<T> *in_out );

// Normalize by RSS (real and complex types only)
template<class REAL, class T> 
bool cuNDA_rss_normalize( cuNDArray<T> *in_out, unsigned int dim );

// Scale (with constant) - real and complex types
template<class A, class X> 
void cuNDA_scale( A a, cuNDArray<X> *x );

// Scale (component-wise) - real and complex types
template<class A, class X> 
bool cuNDA_scale( cuNDArray<A> *a, cuNDArray<X> *x );

// 'axpby' - for real and complex types
template<class A, class B, class XY> 
bool cuNDA_axpby( cuNDArray<A> *a, cuNDArray<XY> *x, cuNDArray<B> *b, cuNDArray<XY> *y );

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

//
// Conversion between vector<unsigned int> and uintd
//

template<unsigned int D> 
std::vector<unsigned int> cuNDA_toVec( typename uintd<D>::Type dims );

template<unsigned int D> 
bool cuNDA_fromVec( std::vector<unsigned int> from, typename uintd<D>::Type &to );
