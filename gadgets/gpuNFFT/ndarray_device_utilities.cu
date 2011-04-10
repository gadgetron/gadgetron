#include "ndarray_device_utilities.hcu"
#include "uintd_operators.hcu"
#include "floatd_operators.hcu"
#include "vector_utilities.hcu"
#include "check_CUDA.h"

#include <cublas.h>
#include <vector>

using namespace std;

// Clear
template<class T> __global__ 
void cuNDA_clear_kernel( T *data, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    T zero = get_zero<T>();
    data[idx] = zero;
  }
}

// Clear
template<class T> __host__
void cuNDA_clear( cuNDArray<T> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make clear image
  cuNDA_clear_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

// Abs
template<class T> __global__ 
void cuNDA_abs_kernel( T *in_out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    T val = in_out[idx];    
    in_out[idx] = abs(val);
  }
}

// Abs
template<class T> __host__ 
void cuNDA_abs( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  // Make modulus image
  cuNDA_abs_kernel<<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

// Reciprocal
template<class T> __global__ 
void cuNDA_reciprocal_kernel( T *data, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    data[idx] = reciprocal(data[idx]);
  }
}

// Reciprocal
template<class T> __host__
void cuNDA_reciprocal( cuNDArray<T> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_reciprocal_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

// Normalize
template<> __host__
void cuNDA_normalize<float, float>( cuNDArray<float> *data, float new_max )
{
  unsigned int number_of_elements = data->get_number_of_elements();

  // Find the maximum value in the array
  int max_idx = cublasIsamax (number_of_elements, data->get_data_ptr(), 1);

  // Copy that value back to host memory
  float max_val;
  cudaMemcpy(&max_val, (data->get_data_ptr()+max_idx-1), sizeof(float), cudaMemcpyDeviceToHost);

  // Scale the array
  cublasSscal( number_of_elements, new_max/max_val, data->get_data_ptr(), 1 );

  CHECK_FOR_CUDA_ERROR();
}

// Normalize
template<> __host__
void cuNDA_normalize<float, cuFloatComplex>( cuNDArray<cuFloatComplex> *data, float new_max )
{
  unsigned int number_of_elements = data->get_number_of_elements();

  // Normalize from magnitudes image
  auto_ptr< cuNDArray<float> > normed_image = cuNDA_norm<float, cuFloatComplex>( data );

  // Find the maximum value in the array
  int max_idx = cublasIsamax( number_of_elements, normed_image->get_data_ptr(), 1 );

  // Copy that value back to host memory
  float max_val;
  cudaMemcpy(&max_val, (normed_image->get_data_ptr()+max_idx-1), sizeof(float), cudaMemcpyDeviceToHost);

  // Scale the array
  cublasCscal( number_of_elements, make_cuFloatComplex(new_max/max_val, 0.0f), data->get_data_ptr(), 1 );

  CHECK_FOR_CUDA_ERROR();
}

// Normalized RSS
template<class REAL, class T> __global__ void
cuNDA_rss_normalize_kernel( T *in_out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

    REAL rss = get_zero<REAL>();
    
    for( unsigned int i=0; i<number_of_batches; i++ )
      rss += norm_sq(in_out[i*stride+in_idx]);

    rss = sqrt(rss);
    rss += get_epsilon<REAL>(); // avoid potential division by zero
    reciprocal(rss);
    
    for( unsigned int i=0; i<number_of_batches; i++ ) {
      T out = in_out[i*stride+in_idx];
      out *= rss;
      in_out[i*stride+in_idx] = out; 
    } 
  }
}

// Normalized RSS
template<class REAL, class T> __host__
bool cuNDA_rss_normalize( cuNDArray<T> *in_out, unsigned int dim )
{
  if( !(in_out->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_rss_normalized:: underdimensioned." << endl; 
    return false;
  }
  
  if( dim > in_out->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_rss_normalized:: dimension out of range." << endl; 
    return false;
  }

  unsigned int number_of_batches = in_out->get_size(dim);
  unsigned int number_of_elements = in_out->get_number_of_elements()/number_of_batches;

  unsigned int stride = 1;
  for( unsigned int i=0; i<dim; i++ )
    stride *= in_out->get_size(i);

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_rss_normalize_kernel<REAL, T><<< gridDim, blockDim >>>( in_out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  return true;
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

// Norm
template<class REAL, class T> __global__ 
void cuNDA_norm_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    T val = in[idx];    
    out[idx] = norm(val);
  }
}

// Norm
template<class REAL, class T> __host__ 
auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<T> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
  
  // Make modulus image
  if( out != 0x0 )
    cuNDA_norm_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<REAL> >(out);
}

// Norm sqaured
template<class REAL, class T> __global__ 
void cuNDA_norm_squared_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    T val = in[idx];    
    out[idx] = norm_sq(val);
  }
}

// Norm squared
template<class REAL, class T> __host__
auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<T> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
  
  // Make norm image
  if( out != 0x0 )
    cuNDA_norm_squared_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<REAL> >(out);
}

// RSS
template<class REAL, class T> __global__ void
cuNDA_rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

    REAL rss = get_zero<REAL>();
    
    for( unsigned int i=0; i<number_of_batches; i++ ) 
      rss += norm_sq(in[i*stride+in_idx]);

    rss = sqrt(rss);

    out[idx] = rss; 
  }
}

// RSS
template<class REAL, class T> __host__ 
auto_ptr< cuNDArray<REAL> > cuNDA_rss( cuNDArray<T> *in, unsigned int dim )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_rss:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<REAL> >(0x0);
  }
  
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_rss:: dimension out of range." << endl; 
    return auto_ptr< cuNDArray<REAL> >(0x0);
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;
  unsigned int stride = 1;

  vector<unsigned int> dims;
  for( unsigned int i=0; i<in->get_number_of_dimensions(); i++ ){
    if( i != dim )
      dims.push_back(in->get_size(i));
    if( i < dim )
      stride *= in->get_size(i);
  }

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(dims);
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if ( out != 0x0 )
    cuNDA_rss_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<REAL> >(out);
}

// Sum
template<class T> __global__ void
cuNDA_sum_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
    
    T val = get_zero<T>();
    
    for( unsigned int i=0; i<number_of_batches; i++ ) 
      val += in[i*stride+in_idx];

    out[idx] = val; 
  }
}

// Sum
template<class T> __host__ 
auto_ptr< cuNDArray<T> > cuNDA_sum( cuNDArray<T> *in, unsigned int dim )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_sum:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<T> >(0x0);
  }
  
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_rss:: dimension out of range." << endl; 
    return auto_ptr< cuNDArray<T> >(0x0);
  }
  
  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;
  unsigned int stride = 1;

  vector<unsigned int> dims;
  for( unsigned int i=0; i<in->get_number_of_dimensions(); i++ ){
    if( i != dim )
      dims.push_back(in->get_size(i));
    if( i < dim )
      stride *= in->get_size(i);
  }

  cuNDArray<T> *out = cuNDArray<T>::allocate(dims);
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if( out != 0x0 )
    cuNDA_sum_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();

  return auto_ptr< cuNDArray<T> >(out);
}

// Build correlation matrix
template<class T> __global__ void
cuNDA_correlation_kernel( T *in, T *corrm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.y;

  if( p < num_elements ){    
    for( unsigned int j=0; j<i; j++){
      corrm[(j*num_batches+i)*num_elements+p] = in[i*num_elements+p]*conj(in[j*num_elements+p]);
      corrm[(i*num_batches+j)*num_elements+p] = conj(corrm[(j*num_batches+i)*num_elements+p]);
    }
    corrm[(i*num_batches+i)*num_elements+p] = in[i*num_elements+p]*conj(in[i*num_elements+p]);
  }
}

// Build correlation matrix
template<class T> __host__ 
auto_ptr< cuNDArray<T> > cuNDA_correlation( cuNDArray<T> *in )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_correlation:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<T> >(0x0);
  }

  unsigned int number_of_batches = in->get_size(in->get_number_of_dimensions()-1);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  vector<unsigned int> dims = in->get_dimensions();
  dims.push_back(number_of_batches);

  cuNDArray<T> *out = cuNDArray<T>::allocate(dims);
  
  int device; cudaGetDevice( &device );
  cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
  unsigned int warp_size = deviceProp.warpSize;

  dim3 blockDim(((512/number_of_batches)/warp_size)*warp_size, number_of_batches);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if( blockDim.x == 0 ){
    cout << endl << "cuNDA_correlation:: correlation dimension exceeds capacity." << endl; 
    return auto_ptr< cuNDArray<T> >(0x0);
  }

  if( out != 0x0 )
    cuNDA_correlation_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<T> >(out);
}

// Crop
template<class UINTd, class T> __global__ void
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

// Crop
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

// Expand and zero fill
template<class UINTd, class T> __global__ void
cuNDA_expand_with_zero_fill_kernel( UINTd matrix_size_in, UINTd matrix_size_out, T *in, T *out, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const UINTd co_out = idx_to_co( idx, matrix_size_out );
    const UINTd offset = (matrix_size_out-matrix_size_in)>>1;
    T _out;
    bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      if( inside )
	_out = in[co_to_idx(co_out-offset, matrix_size_in)+batch*prod(matrix_size_in)];
      else{
	T zero = get_zero<T>();
	_out = zero;
      } 
      out[idx+batch*number_of_elements] = _out;
    }
  }
}

// Expand and zero fill
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

// Zero fill border (rectangular)
template<class UINTd, class T> __global__ void
cuNDA_zero_fill_border_kernel( UINTd matrix_size_in, UINTd matrix_size_out, T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const UINTd co_out = idx_to_co( idx, matrix_size_out );
    const UINTd offset = (matrix_size_out-matrix_size_in)>>1;
    if( weak_less( co_out, offset ) || weak_greater_equal( co_out, matrix_size_in+offset ) ){
      T zero = get_zero<T>();
      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
	image[idx+batch*number_of_elements] = zero;
      }
    }
    else
      ; // do nothing
  }
}

// Zero fill border (rectangular)
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

// Zero fill border (circular)
template<class UINTd, class REALd, class T> __global__ void
cuNDA_zero_fill_border_kernel( REALd radius, UINTd matrix_size_out, T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  
  if( idx < number_of_elements ){
    const UINTd co_out = idx_to_co( idx, matrix_size_out );
    const REALd co_f = abs(uintd_to_reald(co_out) - uintd_to_reald(matrix_size_out>>1));
    if( co_f<radius )
      ; // do nothing
    else{
      T zero = get_zero<T>();
      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
	image[idx+batch*number_of_elements] = zero;
      }
    }    
  }
}

// Zero fill border (circular)
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

// cuNDArray to std::vector
template<class UINTd>  __host__
vector<unsigned int> cuNDA_toVec( UINTd dims )
{
  unsigned int d = sizeof(UINTd)/sizeof(unsigned int);  
  vector<unsigned int> out(d);
  for( unsigned int i=0; i<d; i++ )
    out[i] = ((unsigned int*)&dims)[i];
  return out;
}

// std::vector to cuNDArray
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

// Instanciation

template void cuNDA_clear<float>(cuNDArray<float>*);
template void cuNDA_clear<cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template void cuNDA_abs<float>(cuNDArray<float>*);
template void cuNDA_abs<cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template void cuNDA_reciprocal<float>(cuNDArray<float>*);
template void cuNDA_reciprocal<cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template bool cuNDA_rss_normalize<float, float>(cuNDArray<float>*, unsigned int);
template bool cuNDA_rss_normalize<float, cuFloatComplex>(cuNDArray<cuFloatComplex>*, unsigned int);

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

template auto_ptr< cuNDArray<float> > cuNDA_norm<float, cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float, float>(cuNDArray<float>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float, cuFloatComplex>(cuNDArray<cuFloatComplex>*);

template auto_ptr< cuNDArray<float> > cuNDA_rss<float, float>(cuNDArray<float>*, unsigned int);
template auto_ptr< cuNDArray<float> > cuNDA_rss<float, cuFloatComplex>(cuNDArray<cuFloatComplex>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_sum<float>(cuNDArray<float>*, unsigned int);
template auto_ptr< cuNDArray<cuFloatComplex> > cuNDA_sum<cuFloatComplex>(cuNDArray<cuFloatComplex>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_correlation<float>(cuNDArray<float>*);
template auto_ptr< cuNDArray<cuFloatComplex> > cuNDA_correlation<cuFloatComplex>(cuNDArray<cuFloatComplex>*);

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

