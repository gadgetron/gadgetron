#include "ndarray_device_utilities.hcu"

#include <uintd_operators.hcu>
#include <uintd_utilities.hcu>
#include <floatd_utilities.hcu>
#include "check_CUDA.h"

#include <cublas.h>

#include <vector>

using namespace std;

// Abs (component-wise)
template<class REALd> __global__ 
void cuNDA_abs_kernel( REALd *in_out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    REALd val = in_out[idx];    
    in_out[idx] = abs(val);
  }
}

// Abs (component-wise)
template<class REALd> __host__ 
void cuNDA_abs( cuNDArray<REALd> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  // Make modulus image
  cuNDA_abs_kernel<<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

// Norm
template<class REALd, class REAL> __global__ 
void cuNDA_norm_kernel( REALd *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    REALd val = in[idx];    
    out[idx] = norm(val);
  }
}

// Norm
template<class REALd, class REAL> __host__ 
cuNDArray<REAL>* cuNDA_norm( cuNDArray<REALd> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
  if( out == 0x0 ) return 0x0;
  
  // Make modulus image
  cuNDA_norm_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return out;
}

// Norm sqaured
template<class REALd, class REAL> __global__ 
void cuNDA_norm_squared_kernel( REALd *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    REALd val = in[idx];    
    out[idx] = norm_sq(val);
  }
}

// Norm squared
template<class REALd, class REAL> __host__
cuNDArray<REAL>* cuNDA_norm_squared( cuNDArray<REALd> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
  if( out == 0x0 ) return 0x0;
  
  // Make norm image
  cuNDA_norm_squared_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return out;
}

// Reciprocal
template<class REALd> __global__ 
void cuNDA_reciprocal_kernel( REALd *data, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    data[idx] = reciprocal(data[idx]);
  }
}

// Reciprocal
template<class REALd> __host__
void cuNDA_reciprocal( cuNDArray<REALd> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_reciprocal_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

// Clear
template<class REALd> __global__ 
void cuNDA_clear_kernel( REALd *data, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    REALd zero; get_zero( zero );
    data[idx] = zero;
  }
}

// Clear
template<class REALd> __host__
void cuNDA_clear( cuNDArray<REALd> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make clear image
  cuNDA_clear_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

// Normalize
template<class REAL> __host__
void cuNDA_normalize( cuNDArray<REAL> *data, REAL new_max )
{
  unsigned int number_of_elements = data->get_number_of_elements();

  // Find the maximum value in the array
  int max_idx = cublasIsamax (number_of_elements, data->get_data_ptr(), 1);

  // Copy that value back to host memory
  REAL max_val;
  cudaMemcpy(&max_val, (data->get_data_ptr()+max_idx-1), sizeof(REAL), cudaMemcpyDeviceToHost);

  // Scale the array (2.5 is an "arbitrary" scaling constant)
  cublasSscal( number_of_elements, new_max/max_val, data->get_data_ptr(), 1 );

  CHECK_FOR_CUDA_ERROR();
}

// Scale
template<class A, class X> __global__ 
void cuNDA_scale_kernel( A a, X *x, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    X in = x[idx];
    in *= a;
    x[idx] = in;
  }
}

// Scale 
template<class A, class X> __host__
void cuNDA_scale( A a, cuNDArray<X> *x )
{
  unsigned int number_of_elements = x->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_scale_kernel<<< gridDim, blockDim >>> ( a, x->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

// Scale
template<class A, class X> __global__ 
void cuNDA_scale_kernel( A *a, X *x, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    A in_a = a[idx];
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      X in_x = x[batch*number_of_elements+idx];
      x[batch*number_of_elements+idx] = in_a*in_x;
    }
  }
}

// Scale 
template<class A, class X> __host__
bool cuNDA_scale( cuNDArray<A> *a, cuNDArray<X> *x )
{
  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch, cannot scale" << endl;
    return false;
  }
  
  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  // Invoke kernel
  cuNDA_scale_kernel<<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), num_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// 'axpy'
template<class A, class XY> __global__ 
void cuNDA_axpy_kernel( A a, XY *x, XY *y, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    XY in_x = x[idx];
    XY in_y = y[idx];
    in_y += (a*in_x);
    y[idx] = in_y;
  }
}

// 'axpy' 
template<class A, class XY> __host__
bool cuNDA_axpy( A a, cuNDArray<XY> *x, cuNDArray<XY> *y )
{
  if( x->get_number_of_elements() != y->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpy'" << endl;
    return false;
  }

  unsigned int number_of_elements = y->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_axpy_kernel<<< gridDim, blockDim >>> ( a, x->get_data_ptr(), y->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// '.axpby'
template<class A, class B, class XY> __global__ 
void cuNDA_axpby_kernel( A *a, XY *x, B *b, XY *y, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    A in_a = a[idx];
    B in_b = b[idx];
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      unsigned int iidx = batch*number_of_elements + idx;
      XY in_x = x[iidx];
      XY in_y = y[iidx];
      in_y *= in_b;
      in_y += (in_a*in_x);
      y[iidx] = in_y;
    }
  }
}

// '.axpby' 
template<class A, class B, class XY> __host__
bool cuNDA_axpby( cuNDArray<A> *a, cuNDArray<XY> *x, cuNDArray<B> *b, cuNDArray<XY> *y )
{

  if( x->get_number_of_elements() != y->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpby'" << endl;
    return false;
  }

  if( a->get_number_of_elements() != b->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpby'" << endl;
    return false;
  }

  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpby'" << endl;
    return false;
  }
  
  unsigned int number_of_batches = x->get_number_of_elements() / a->get_number_of_elements();
  unsigned int number_of_elements = a->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_axpby_kernel<<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), b->get_data_ptr(), y->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// Crop
template< class UINTd, class T > __global__ void
cuNDA_crop_kernel( UINTd offset, UINTd matrix_size_in, UINTd matrix_size_out, T *in, T *out, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const UINTd co = idx_to_co( idx, matrix_size_out );
    const UINTd co_os = offset + co;
    const unsigned int source_idx = co_to_idx(co_os, matrix_size_in);
    const unsigned int source_elements = prod(matrix_size_in);
    for( unsigned int image=0; image<number_of_batches; image++ )
      out[image*number_of_elements+idx] = in[image*source_elements+source_idx];
  }
}

template<class UINTd, class T> __host__
bool cuNDA_crop( UINTd offset, cuNDArray<T> *in, cuNDArray<T> *out )
{ 
  unsigned int d = sizeof(UINTd)/sizeof(unsigned int);

  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "image dimensions mismatch, cannot crop" << endl;
    return false;
  }

  if( !(in->get_number_of_dimensions() == d || in->get_number_of_dimensions() == d+1) ){
    cout << endl << "image dimensions mismatch, cannot crop" << endl;
    return false;
  }

  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == d ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  UINTd matrix_size_in; cuNDA_fromVec( in->get_dimensions(), matrix_size_in );
  UINTd matrix_size_out; cuNDA_fromVec( out->get_dimensions(), matrix_size_out );
  
  if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
    cout << endl << "cropping size mismatch, cannot crop" << endl;
    return false;
  }
    
  unsigned int number_of_elements = prod(matrix_size_out);

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_crop_kernel<<< gridDim, blockDim >>> ( offset, matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}

// Expand_With_Zero_Fill
template< class UINTd, class T > __global__ void
cuNDA_expand_with_zero_fill_kernel( UINTd matrix_size_in, UINTd matrix_size_out, T *in, T *out, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const UINTd co_out = idx_to_co( idx, matrix_size_out );
    const UINTd offset = (matrix_size_out-matrix_size_in)>>1;
    T _out;
    bool inside = weak_greater_equal( co_out, offset ) && weak_less( co_out, matrix_size_in+offset );
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      if( inside )
	_out = in[co_to_idx(co_out-offset, matrix_size_in)+batch*prod(matrix_size_in)];
      else{
	T zero; get_zero(zero);
	_out = zero;
      } 
      out[idx+batch*number_of_elements] = _out;
    }
  }
}

template<class UINTd, class T> __host__
bool cuNDA_expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out )
{ 
  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "Image dimensions mismatch, cannot expand" << endl;
    return false;
  }
  
  UINTd matrix_size_in; cuNDA_fromVec( in->get_dimensions(), matrix_size_in );
  UINTd matrix_size_out; cuNDA_fromVec( out->get_dimensions(), matrix_size_out );
  
  if( weak_greater(matrix_size_in,matrix_size_out) ){
    cout << endl << "Size mismatch, cannot expand" << endl;
    return false;
  }
    
  unsigned int d = sizeof(UINTd)/sizeof(unsigned int);

  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == d ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  unsigned int number_of_elements = prod(matrix_size_out);
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  // Invoke kernel
  cuNDA_expand_with_zero_fill_kernel<<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}

// Zero_Fill_Border
template< class UINTd, class T > __global__ void
cuNDA_zero_fill_border_kernel( UINTd matrix_size_in, UINTd matrix_size_out, T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const UINTd co_out = idx_to_co( idx, matrix_size_out );
    const UINTd offset = (matrix_size_out-matrix_size_in)>>1;
    if( weak_less( co_out, offset ) || weak_greater_equal( co_out, matrix_size_in+offset ) ){
      T zero; get_zero(zero);
      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
	image[idx+batch*number_of_elements] = zero;
      }
    }
    else
      ; // do nothing
  }
}

template<class UINTd, class T> __host__
bool cuNDA_zero_fill_border( UINTd matrix_size_in, cuNDArray<T> *out )
{ 
  unsigned int d = sizeof(UINTd)/sizeof(unsigned int);
  
  UINTd matrix_size_out; cuNDA_fromVec( out->get_dimensions(), matrix_size_out );
  
  if( weak_greater(matrix_size_in, matrix_size_out) ){
    cout << endl << "Size mismatch, cannot zero fill" << endl;
    return false;
  }
    
  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == d ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  unsigned int number_of_elements = prod(matrix_size_out);
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  // Invoke kernel
  cuNDA_zero_fill_border_kernel<<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, out->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}

// Zero_Fill_Border
template< class UINTd, class REALd, class T > __global__ void
cuNDA_zero_fill_border_kernel( REALd radius, UINTd matrix_size_out, T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  
  if( idx < number_of_elements ){
    const UINTd co_out = idx_to_co( idx, matrix_size_out );
    const REALd co_f = abs(uintd_to_reald(co_out) - uintd_to_reald(matrix_size_out));
    if( weak_less( co_f, radius ) )
      ; // do nothing
    else{
      T zero; get_zero(zero);
      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
	image[idx+batch*number_of_elements] = zero;
      }
    }    
  }
}

template<class UINTd, class REALd, class T> __host__
bool cuNDA_zero_fill_border( REALd radius, cuNDArray<T> *out )
{
  unsigned int d = sizeof(UINTd)/sizeof(unsigned int);
  
  if( out->get_number_of_dimensions() != d ){
    cout << endl << "Image dimensions mismatch, cannot zero fill" << endl;
    return false;
  }
  
  UINTd matrix_size_out; cuNDA_fromVec( out->get_dimensions(), matrix_size_out );
  
  if( weak_greater(radius, huintd_to_reald(matrix_size_out)) ){
    cout << endl << "Size mismatch, cannot zero fill" << endl;
    return false;
  }
  
  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == d ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  unsigned int number_of_elements = prod(matrix_size_out);
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  // Invoke kernel
  cuNDA_zero_fill_border_kernel<<< gridDim, blockDim >>> ( radius, matrix_size_out, out->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return true;
}

template<class UINTd>  __host__
vector<unsigned int> cuNDA_toVec( UINTd dims )
{
  unsigned int d = sizeof(UINTd)/sizeof(unsigned int);  
  vector<unsigned int> out(d);
  for( unsigned int i=0; i<d; i++ )
    out[i] = ((unsigned int*)&dims)[i];
  return out;
}

template<class UINTd>  __host__
bool cuNDA_fromVec( vector<unsigned int> from, UINTd &to )
{
  unsigned int d = sizeof(UINTd)/sizeof(unsigned int); 
  if( from.size() < d ){
    cout << "Cannot convert vector to UINTd" << endl;
    return false;
  }
  
  vector<unsigned int>::iterator it = from.begin();
  for( unsigned int i=0; i<d; i++ )
    ((unsigned int*)&to)[i] = *it;

  return true;
}


/*
template< class T > __global__ void 
add_images_kernel( unsigned int num_elements, T *target, T *source1, T *source2  )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < num_elements )
    target[idx] = __image_add( source1[idx], source2[idx] );
}

template< class T > __host__ void 
add_images( unsigned int num_elements, T *targetDevPtr, T *source1DevPtr,  T *source2DevPtr  )
{

  // Find dimensions of grid/blocks.

  cudaDeviceProp deviceProp;  
  cudaGetDeviceProperties( &deviceProp, _convolution_device );

  dim3 blockDim( deviceProp.maxThreadsPerBlock, 1, 1 );
  dim3 gridDim( (unsigned int) ceil((double)num_elements/blockDim.x), 1, 1 );

  // Invoke kernel
  add_images_kernel<T><<< gridDim, blockDim >>>( num_elements, targetDevPtr, source1DevPtr, source2DevPtr );

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    printf("\nCuda error detected in 'add_images_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
    exit(1);
  }
}

template< class T > __global__ void 
add_images_kernel( unsigned int num_elements, T *target, T *source, unsigned int num_images  )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < num_elements ){

    T out;

    __zero_element( &out );

    for( unsigned int i=0; i<num_images; i++ )
      out = __image_add( out, source[i*num_elements+idx] );

    target[idx] = out;
  }
}

template< class T > void 
add_images( unsigned int num_elements, T *targetDevPtr, T *sourceDevPtr, unsigned int num_source_images )
{

  // Find dimensions of grid/blocks.

  cudaDeviceProp deviceProp;  
  cudaGetDeviceProperties( &deviceProp, _convolution_device );

  dim3 blockDim( deviceProp.maxThreadsPerBlock, 1, 1 );
  dim3 gridDim( (unsigned int) ceil((double)num_elements/blockDim.x), 1, 1 );

  // Invoke kernel
  add_images_kernel<T><<< gridDim, blockDim >>>( num_elements, targetDevPtr, sourceDevPtr, num_source_images );

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    printf("\nCuda error detected in 'add_images_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
    exit(1);
  }
}
*/

// Instanciation

template void cuNDA_clear<float>(cuNDArray<float>*);
template void cuNDA_clear<cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template void cuNDA_abs<float>(cuNDArray<float>*);
template void cuNDA_abs<cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template cuNDArray<float>* cuNDA_norm<cuFloatComplex, float>(cuNDArray<cuFloatComplex>*);

template cuNDArray<float>* cuNDA_norm_squared<cuFloatComplex, float>(cuNDArray<cuFloatComplex>*);

template void cuNDA_reciprocal<float>(cuNDArray<float>*);
template void cuNDA_reciprocal<cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template void cuNDA_normalize<float>( cuNDArray<float>*, float);

template void cuNDA_scale<float, float>(float, cuNDArray<float>*);
template void cuNDA_scale<float, cuFloatComplex>(float, cuNDArray<cuFloatComplex>*);
template void cuNDA_scale<cuFloatComplex, cuFloatComplex>(cuFloatComplex, cuNDArray<cuFloatComplex>*);

template bool cuNDA_scale<float, float>(cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_scale<float, cuFloatComplex>(cuNDArray<float>*, cuNDArray<cuFloatComplex>*);
template bool cuNDA_scale<cuFloatComplex, cuFloatComplex>(cuNDArray<cuFloatComplex>*, cuNDArray<cuFloatComplex>*);

template bool cuNDA_axpy<float, float>( float, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_axpy<float, cuFloatComplex>( float, cuNDArray<cuFloatComplex>*, cuNDArray<cuFloatComplex>*);

template bool cuNDA_axpby<float, float, float>( cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_axpby<float, float, cuFloatComplex>( cuNDArray<float>*, cuNDArray<cuFloatComplex>*, cuNDArray<float>*, cuNDArray<cuFloatComplex>*);

template bool cuNDA_crop<uint2, float>(uint2, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<uint2, cuFloatComplex>(uint2, cuNDArray<cuFloatComplex>*, cuNDArray<cuFloatComplex>*);

template bool cuNDA_expand_with_zero_fill<uint2, float>(cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<uint2, cuFloatComplex>(cuNDArray<cuFloatComplex>*, cuNDArray<cuFloatComplex>*);

template bool cuNDA_zero_fill_border<uint2, float>(uint2, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<uint2, cuFloatComplex>(uint2, cuNDArray<cuFloatComplex>*);

template bool cuNDA_zero_fill_border<uint2, float2, float>(float2, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<uint2, float2, cuFloatComplex>(float2, cuNDArray<cuFloatComplex>*);

template std::vector<unsigned int> cuNDA_toVec<uint2>(uint2);

template bool cuNDA_fromVec<uint2>( std::vector<unsigned int>, uint2&);

