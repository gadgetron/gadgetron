#include "ndarray_vector_td_utilities.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"

#include <vector>

using namespace std;

// Sum
template<class T> __global__ void
cuNDA_sum_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
 
    T val = in[in_idx];
 
    for( unsigned int i=1; i<number_of_batches; i++ ) 
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
    cout << endl << "cuNDA_sum:: dimension out of range." << endl; 
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
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if( out != 0x0 )
    cuNDA_sum_kernel<T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  return auto_ptr< cuNDArray<T> >(out);
}

// Norm
template<class REAL, unsigned int D> __global__ 
void cuNDA_norm_kernel( typename reald<REAL,D>::Type *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    typename reald<REAL,D>::Type val = in[idx]; 
    out[idx] = norm<REAL,D>(val);
  }
}

// Norm
template<class REAL> __host__ 
auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<REAL> *in )
{
  return cuNDA_norm<REAL,1>((cuNDArray<typename reald<REAL,1>::Type>*) in);
}

// Norm
template<class REAL, unsigned int D> __host__ 
auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<typename reald<REAL,D>::Type> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
 
  // Make modulus image
  if( out != 0x0 )
    cuNDA_norm_kernel<REAL,D><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<REAL> >(out);
}

// Norm sqaured
template<class REAL, unsigned int D> __global__ 
void cuNDA_norm_squared_kernel( typename reald<REAL,D>::Type *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    typename reald<REAL,D>::Type val = in[idx]; 
    out[idx] = norm_squared<REAL,D>(val);
  }
}

// Norm squared
template<class REAL> __host__
auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<REAL> *in )
{
  return cuNDA_norm_squared<REAL,1>((cuNDArray<typename reald<REAL,1>::Type>*) in);
}

// Norm squared
template<class REAL, unsigned int D> __host__
auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<typename reald<REAL,D>::Type> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
 
  // Make norm image
  if( out != 0x0 )
    cuNDA_norm_squared_kernel<REAL,D><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
 
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
      rss += norm_squared<REAL>(in[i*stride+in_idx]);

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
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if ( out != 0x0 )
    cuNDA_rss_kernel<REAL,T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<REAL> >(out);
}

// Build correlation matrix
template<class REAL, class T> __global__ void
cuNDA_correlation_kernel( T *in, T *corrm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.y;

  if( p < num_elements ){
    for( unsigned int j=0; j<i; j++){
      T tmp = mul<T,T>(in[i*num_elements+p], conj<REAL>(in[j*num_elements+p]));
      corrm[(j*num_batches+i)*num_elements+p] = tmp;
      corrm[(i*num_batches+j)*num_elements+p] = conj<REAL>(tmp);
    }
    T tmp = in[i*num_elements+p];
    corrm[(i*num_batches+i)*num_elements+p] = mul<T,T>(tmp,conj<REAL>(tmp));
  }
}

// Build correlation matrix
template<class REAL> __host__ 
auto_ptr< cuNDArray<typename complext<REAL>::Type> > cuNDA_correlation( cuNDArray<typename complext<REAL>::Type> *in )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_correlation:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
  }

  unsigned int number_of_batches = in->get_size(in->get_number_of_dimensions()-1);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  vector<unsigned int> dims = in->get_dimensions();
  dims.push_back(number_of_batches);

  cuNDArray<typename complext<REAL>::Type> *out = cuNDArray<typename complext<REAL>::Type>::allocate(dims);
 
  int device; cudaGetDevice( &device );
  cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
  unsigned int warp_size = deviceProp.warpSize;

  dim3 blockDim(((deviceProp.maxThreadsPerBlock/number_of_batches)/warp_size)*warp_size, number_of_batches);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if( blockDim.x == 0 ){
    cout << endl << "cuNDA_correlation:: correlation dimension exceeds capacity." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
  }

  if( out != 0x0 )
    cuNDA_correlation_kernel<REAL><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(out);
}

// Clear
template<class T> __global__ 
void cuNDA_clear_kernel( T *in_out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T zero = get_zero<T>();
    in_out[idx] = zero;
  }
}

// Clear
template<class T> __host__
void cuNDA_clear( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make clear image
  cuNDA_clear_kernel<<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
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

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
 
  // Make modulus image
  cuNDA_abs_kernel<<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
}

// Reciprocal
template<class T> __global__ 
void cuNDA_reciprocal_kernel( T *in_out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = reciprocal<T>(in_out[idx]);
  }
}

// Reciprocal
template<class T> __host__
void cuNDA_reciprocal( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_reciprocal_kernel<<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
}

// Normalize (float)
__host__
void cuNDA_normalize( cuNDArray<float> *data, float new_max, cublasHandle_t handle )
{
  unsigned int number_of_elements = data->get_number_of_elements();

  // Find the maximum value in the array
  int max_idx;
  cublasIsamax( handle, number_of_elements, data->get_data_ptr(), 1, &max_idx );

  // Copy that value back to host memory
  float max_val;
  cudaMemcpy(&max_val, (data->get_data_ptr()+max_idx-1), sizeof(float), cudaMemcpyDeviceToHost);

  // Scale the array
  float scale = new_max/max_val;
  cublasSscal( handle, number_of_elements, &scale, data->get_data_ptr(), 1 );

  CHECK_FOR_CUDA_ERROR();
}

// Normalize (double)
__host__
void cuNDA_normalize( cuNDArray<double> *data, double new_max, cublasHandle_t handle )
{
  unsigned int number_of_elements = data->get_number_of_elements();

  // Find the maximum value in the array
  int max_idx;
  cublasIdamax( handle, number_of_elements, data->get_data_ptr(), 1, &max_idx );

  // Copy that value back to host memory
  double max_val;
  cudaMemcpy(&max_val, (data->get_data_ptr()+max_idx-1), sizeof(double), cudaMemcpyDeviceToHost);

  // Scale the array
  double scale = new_max/max_val;
  cublasDscal( handle, number_of_elements, &scale, data->get_data_ptr(), 1 );

  CHECK_FOR_CUDA_ERROR();
}

// Normalize
template<class REAL> __host__
void cuNDA_normalize( cuNDArray<REAL> *data, REAL new_max )
{
  cuNDA_normalize( data, new_max );
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
      rss += norm_squared<REAL>(in_out[i*stride+in_idx]);
 
    rss = sqrt(rss);
    rss += get_epsilon<REAL>(); // avoid potential division by zero
    rss = reciprocal(rss);
 
    for( unsigned int i=0; i<number_of_batches; i++ ) {
      T out = in_out[i*stride+in_idx];
      out *= rss; // complex-scalar multiplication (element-wise operator)
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

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_rss_normalize_kernel<REAL,T><<< gridDim, blockDim >>>
    ( in_out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// Scale
template<class A, class X> __global__ 
void cuNDA_scale1_kernel( A a, X *x, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    X in = x[idx];
    in = mul<A,X>(a,in);
    x[idx] = in;
  }
}

// Scale 
template<class A, class X> __host__
void cuNDA_scale( A a, cuNDArray<X> *x )
{
  unsigned int number_of_elements = x->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_scale1_kernel<A,X><<< gridDim, blockDim >>> ( a, x->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
}

// Scale
template<class A, class X> __global__ 
void cuNDA_scale2_kernel( A *a, X *x, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    A in_a = a[idx];
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      X in_x = x[batch*number_of_elements+idx];
      x[batch*number_of_elements+idx] = mul<A,X>(in_a,in_x);
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
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
 
  // Invoke kernel
  cuNDA_scale2_kernel<A,X><<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), num_batches, number_of_elements );
 
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
    in_y += mul<A,XY>(a,in_x);
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

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_axpy_kernel<<< gridDim, blockDim >>> ( a, x->get_data_ptr(), y->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// 'axpby'
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
      in_y = mul<B,XY>(in_b,in_y);
      in_y += mul<A,XY>(in_a,in_x);
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

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_axpby_kernel<<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), b->get_data_ptr(), y->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// Crop
template<class T, unsigned int D> __global__ void
cuNDA_crop_kernel( typename uintd<D>::Type offset, typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
		   T *in, T *out, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const typename uintd<D>::Type co = idx_to_co<D>( idx, matrix_size_out );
    const typename uintd<D>::Type co_os = offset + co;
    const unsigned int source_idx = co_to_idx<D>(co_os, matrix_size_in);
    const unsigned int source_elements = prod(matrix_size_in);
    for( unsigned int image=0; image<number_of_batches; image++ )
      out[image*number_of_elements+idx] = in[image*source_elements+source_idx];
  }
}

// Crop
template<class T, unsigned int D> __host__
bool cuNDA_crop( typename uintd<D>::Type offset, cuNDArray<T> *in, cuNDArray<T> *out )
{ 
  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "image dimensions mismatch, cannot crop" << endl;
    return false;
  }

  if( !(in->get_number_of_dimensions() == D || in->get_number_of_dimensions() == D+1) ){
    cout << endl << "image dimensions mismatch, cannot crop" << endl;
    return false;
  }

  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == D ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  typename uintd<D>::Type matrix_size_in; cuNDA_fromVec<D>( in->get_dimensions(), matrix_size_in );
  typename uintd<D>::Type matrix_size_out; cuNDA_fromVec<D>( out->get_dimensions(), matrix_size_out );
 
  if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
    cout << endl << "cropping size mismatch, cannot crop" << endl;
    return false;
  }
 
  unsigned int number_of_elements = prod(matrix_size_out);

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_crop_kernel<T,D><<< gridDim, blockDim >>> ( offset, matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  return true;
}

// Expand and zero fill
template<class T, unsigned int D> __global__ void
cuNDA_expand_with_zero_fill_kernel( typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
				    T *in, T *out, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
    const typename uintd<D>::Type offset = (matrix_size_out-matrix_size_in)>>1;
    T _out;
    bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      if( inside )
	_out = in[co_to_idx<D>(co_out-offset, matrix_size_in)+batch*prod(matrix_size_in)];
      else{
	T zero = get_zero<T>();
	_out = zero;
      } 
      out[idx+batch*number_of_elements] = _out;
    }
  }
}

// Expand and zero fill
template<class T, unsigned int D> __host__
bool cuNDA_expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out )
{ 
  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "Image dimensions mismatch, cannot expand" << endl;
    return false;
  }
 
  typename uintd<D>::Type matrix_size_in; cuNDA_fromVec<D>( in->get_dimensions(), matrix_size_in );
  typename uintd<D>::Type matrix_size_out; cuNDA_fromVec<D>( out->get_dimensions(), matrix_size_out );
 
  if( weak_greater(matrix_size_in,matrix_size_out) ){
    cout << endl << "Size mismatch, cannot expand" << endl;
    return false;
  }
 
  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == D ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  unsigned int number_of_elements = prod(matrix_size_out);
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
 
  // Invoke kernel
  cuNDA_expand_with_zero_fill_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  return true;
}

// Zero fill border (rectangular)
template<class T, unsigned int D> __global__ void
cuNDA_zero_fill_border_kernel( typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
			       T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
    const typename uintd<D>::Type offset = (matrix_size_out-matrix_size_in)>>1;
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
template<class T, unsigned int D> __host__
bool cuNDA_zero_fill_border( typename uintd<D>::Type matrix_size_in, cuNDArray<T> *out )
{ 
  typename uintd<D>::Type matrix_size_out; cuNDA_fromVec<D>( out->get_dimensions(), matrix_size_out );
 
  if( weak_greater(matrix_size_in, matrix_size_out) ){
    cout << endl << "Size mismatch, cannot zero fill" << endl;
    return false;
  }
 
  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == D ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  unsigned int number_of_elements = prod(matrix_size_out);
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
 
  // Invoke kernel
  cuNDA_zero_fill_border_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  return true;
}

// Zero fill border (circular)
template<class REAL, class T, unsigned int D> __global__ void
cuNDA_zero_fill_border_kernel( typename reald<REAL,D>::Type radius, typename uintd<D>::Type matrix_size_out, 
			       T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
 
  if( idx < number_of_elements ){
    
    typename reald<REAL,D>::Type half_matrix_size_out_real; to_reald( half_matrix_size_out_real, matrix_size_out>>1 );

    const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
    typename reald<REAL,D>::Type co_out_real; to_reald( co_out_real, co_out );
    
    typename reald<REAL,D>::Type co_f = abs( co_out_real - half_matrix_size_out_real );
    
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
template<class REAL, class T, unsigned int D> __host__
bool cuNDA_zero_fill_border( typename reald<REAL,D>::Type radius, cuNDArray<T> *out )
{
  if( out->get_number_of_dimensions() != D ){
    cout << endl << "Image dimensions mismatch, cannot zero fill" << endl;
    return false;
  }
 
  typename uintd<D>::Type matrix_size_out; cuNDA_fromVec<D>( out->get_dimensions(), matrix_size_out );
  typename reald<REAL,D>::Type matrix_size_out_real; to_reald( matrix_size_out_real, matrix_size_out );

  if( weak_greater(radius, matrix_size_out_real) ){
    cout << endl << "Size mismatch, cannot zero fill" << endl;
    return false;
  }
 
  unsigned int number_of_batches = 
    (out->get_number_of_dimensions() == D ) ? 1 : out->get_size(out->get_number_of_dimensions()-1);

  unsigned int number_of_elements = prod(matrix_size_out);
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
 
  // Invoke kernel
  cuNDA_zero_fill_border_kernel<REAL,T,D><<< gridDim, blockDim >>> ( radius, matrix_size_out, out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return true;
}

template<class T> T _dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, cublasHandle_t handle );

template<> typename float_complext::Type
_dot( cuNDArray<typename float_complext::Type>* arr1, cuNDArray<typename float_complext::Type>* arr2, cublasHandle_t handle )
{
  typename float_complext::Type ret = get_zero<float_complext::Type>();
  
  if (cublasCdotc( handle, arr1->get_number_of_elements(),
		   (const cuComplex*) arr1->get_data_ptr(), 1, 
		   (const cuComplex*) arr2->get_data_ptr(), 1,
		   (cuComplex*) &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculating using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> float 
_dot( cuNDArray<float>* arr1, cuNDArray<float>* arr2, cublasHandle_t handle )
{
  float ret = get_zero<float>();
  if( cublasSdot(handle, arr1->get_number_of_elements(),
		 arr1->get_data_ptr(), 1, 
		 arr2->get_data_ptr(), 1,
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculating using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> typename double_complext::Type
_dot( cuNDArray<typename double_complext::Type>* arr1, cuNDArray<typename double_complext::Type>* arr2, cublasHandle_t handle )
{
  typename double_complext::Type ret = get_zero<double_complext::Type>();
  
  if( cublasZdotc(handle, arr1->get_number_of_elements(),
		  (const cuDoubleComplex*) arr1->get_data_ptr(), 1, 
		  (const cuDoubleComplex*) arr2->get_data_ptr(), 1,
		  (cuDoubleComplex*) &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculating using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> double 
_dot( cuNDArray<double>* arr1, cuNDArray<double>* arr2, cublasHandle_t handle )
{
  double ret = get_zero<double>();
  if( cublasDdot(handle, arr1->get_number_of_elements(),
		 arr1->get_data_ptr(), 1, 
		 arr2->get_data_ptr(), 1,
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculating using cublas failed" << std::endl;
    }
  
  return ret;
}

template<class T> T
cuNDA_dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, cublasHandle_t handle )
{
  if (arr1->get_number_of_elements() != arr2->get_number_of_elements()) {
    cout << "cuNDA_dot: array dimensions mismatch" << std::endl;
    return get_zero<T>();
  }

  return _dot<T>( arr1, arr2, handle );  
}

template<class T> bool 
_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y, cublasHandle_t handle );

template<> bool 
_axpy( float_complext::Type a, cuNDArray<float_complext::Type>* x, cuNDArray<float_complext::Type>* y, cublasHandle_t handle )
{
  if( cublasCaxpy( handle, x->get_number_of_elements(), (cuComplex*) &a,
		   (cuComplex*) x->get_data_ptr(), 1, 
		   (cuComplex*) y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: axpy calculating using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<> bool
_axpy( float a, cuNDArray<float>* x, cuNDArray<float>* y, cublasHandle_t handle )
{
  if( cublasSaxpy( handle, x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1, 
		   y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: axpy calculating using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<> bool 
_axpy( double_complext::Type a, cuNDArray<double_complext::Type>* x, cuNDArray<double_complext::Type>* y, cublasHandle_t handle )
{
  if( cublasZaxpy( handle, x->get_number_of_elements(), (cuDoubleComplex*) &a,
		   (cuDoubleComplex*) x->get_data_ptr(), 1, 
		   (cuDoubleComplex*) y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: axpy calculating using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<> bool
_axpy( double a, cuNDArray<double>* x, cuNDArray<double>* y, cublasHandle_t handle )
{
  if( cublasDaxpy( handle, x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1, 
		   y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: axpy calculating using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<class T> bool
cuNDA_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y, cublasHandle_t handle )
{
  if (x->get_number_of_elements() != y->get_number_of_elements()) {
    cout << "cuNDA_dot: axpy array dimensions mismatch" << std::endl;
    return false;
  }
  
  return _axpy<T>( a, x, y, handle );
}

template<class T> bool 
_scal( T a, cuNDArray<T>* x, cublasHandle_t handle );

template<> bool
_scal( float_complext::Type a, cuNDArray<float_complext::Type>* x, cublasHandle_t handle) 
{
  if( cublasCscal( handle, x->get_number_of_elements(), (cuComplex*) &a,
		   (cuComplex*) x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "ccuNDA_scal: calculating using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<> bool
_scal( float a, cuNDArray<float>* x, cublasHandle_t handle ) 
{
  if( cublasSscal( handle, x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "ccuNDA_scal: calculating using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<> bool
_scal( double_complext::Type a, cuNDArray<double_complext::Type>* x, cublasHandle_t handle) 
{
  if( cublasZscal( handle, x->get_number_of_elements(), (cuDoubleComplex*) &a,
		   (cuDoubleComplex*) x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "ccuNDA_scal: calculating using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<> bool
_scal( double a, cuNDArray<double>* x, cublasHandle_t handle ) 
{
  if( cublasDscal( handle, x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "ccuNDA_scal: calculating using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<class T> bool
cuNDA_scal( T a, cuNDArray<T>* x, cublasHandle_t handle )
{
  return _scal<T>( a, x, handle );
}

// cuNDArray to std::vector
template<unsigned int D> __host__
vector<unsigned int> cuNDA_toVec( typename uintd<D>::Type dims )
{
  vector<unsigned int> out(D);
  for( unsigned int i=0; i<D; i++ )
    out[i] = dims.vec[i];
  return out;
}

// std::vector to cuNDArray
template<unsigned int D> __host__
bool cuNDA_fromVec( vector<unsigned int> from, typename uintd<D>::Type &to )
{
  if( from.size() < D ){
    cout << "Cannot convert vector to Typename Uintd" << endl;
    return false;
  }
 
  vector<unsigned int>::iterator it = from.begin();
  for( unsigned int i=0; i<D; i++ ){
    to.vec[i] = *it; it++;
  }
 
  return true;
}

//
// Instantiation
//

template std::vector<unsigned int> cuNDA_toVec<1>( typename uintd<1>::Type );
template std::vector<unsigned int> cuNDA_toVec<2>( typename uintd<2>::Type );
template std::vector<unsigned int> cuNDA_toVec<3>( typename uintd<3>::Type );
template std::vector<unsigned int> cuNDA_toVec<4>( typename uintd<4>::Type );

template bool cuNDA_fromVec<1>( std::vector<unsigned int>, typename uintd<1>::Type& );
template bool cuNDA_fromVec<2>( std::vector<unsigned int>, typename uintd<2>::Type& );
template bool cuNDA_fromVec<3>( std::vector<unsigned int>, typename uintd<3>::Type& );
template bool cuNDA_fromVec<4>( std::vector<unsigned int>, typename uintd<4>::Type& );

template auto_ptr< cuNDArray<int> > cuNDA_sum<int>( cuNDArray<int>*, unsigned int);
template auto_ptr< cuNDArray<intd<1>::Type> > cuNDA_sum<intd<1>::Type >( cuNDArray<intd<1>::Type >*, unsigned int );
template auto_ptr< cuNDArray<intd<2>::Type> > cuNDA_sum<intd<2>::Type >( cuNDArray<intd<2>::Type >*, unsigned int );
template auto_ptr< cuNDArray<intd<3>::Type> > cuNDA_sum<intd<3>::Type >( cuNDArray<intd<3>::Type >*, unsigned int );
template auto_ptr< cuNDArray<intd<4>::Type> > cuNDA_sum<intd<4>::Type >( cuNDArray<intd<4>::Type >*, unsigned int );

template auto_ptr< cuNDArray<unsigned int> > cuNDA_sum<unsigned int>( cuNDArray<unsigned int>*, unsigned int);
template auto_ptr< cuNDArray<uintd<1>::Type> > cuNDA_sum<uintd<1>::Type>( cuNDArray<uintd<1>::Type>*, unsigned int );
template auto_ptr< cuNDArray<uintd<2>::Type> > cuNDA_sum<uintd<2>::Type>( cuNDArray<uintd<2>::Type>*, unsigned int );
template auto_ptr< cuNDArray<uintd<3>::Type> > cuNDA_sum<uintd<3>::Type>( cuNDArray<uintd<3>::Type>*, unsigned int );
template auto_ptr< cuNDArray<uintd<4>::Type> > cuNDA_sum<uintd<4>::Type>( cuNDArray<uintd<4>::Type>*, unsigned int );

// Instanciation -- single precision

template auto_ptr< cuNDArray<float> > cuNDA_sum<float>( cuNDArray<float>*, unsigned int);
template auto_ptr< cuNDArray<floatd<1>::Type> > cuNDA_sum<floatd<1>::Type>( cuNDArray<floatd<1>::Type>*, unsigned int );
template auto_ptr< cuNDArray<floatd<2>::Type> > cuNDA_sum<floatd<2>::Type>( cuNDArray<floatd<2>::Type>*, unsigned int );
template auto_ptr< cuNDArray<floatd<3>::Type> > cuNDA_sum<floatd<3>::Type>( cuNDArray<floatd<3>::Type>*, unsigned int );
template auto_ptr< cuNDArray<floatd<4>::Type> > cuNDA_sum<floatd<4>::Type>( cuNDArray<floatd<4>::Type>*, unsigned int );

template auto_ptr< cuNDArray<float> > cuNDA_norm<float>( cuNDArray<float>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm<float,1>( cuNDArray<floatd<1>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm<float,2>( cuNDArray<floatd<2>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm<float,3>( cuNDArray<floatd<3>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm<float,4>( cuNDArray<floatd<4>::Type>*);

template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float>( cuNDArray<float>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,1>( cuNDArray<floatd<1>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,2>( cuNDArray<floatd<2>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,3>( cuNDArray<floatd<3>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,4>( cuNDArray<floatd<4>::Type>*);

template auto_ptr< cuNDArray<float> > cuNDA_rss<float, float>( cuNDArray<float>*, unsigned int);
template auto_ptr< cuNDArray<float> > cuNDA_rss<float, typename complext<float>::Type>( cuNDArray<typename complext<float>::Type>*, unsigned int);

template auto_ptr< cuNDArray<typename complext<float>::Type> > cuNDA_correlation<float>( cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_crop<float,1>( typename uintd<1>::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,1>( typename uintd<1>::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,1>( typename uintd<1>::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,1>( typename uintd<1>::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,1>( typename uintd<1>::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_crop<float,2>( typename uintd<2>::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,2>( typename uintd<2>::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,2>( typename uintd<2>::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,2>( typename uintd<2>::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,2>( typename uintd<2>::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_crop<float,3>( typename uintd<3>::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,3>( typename uintd<3>::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,3>( typename uintd<3>::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,3>( typename uintd<3>::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,3>( typename uintd<3>::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_crop<float,4>( typename uintd<4>::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,4>( typename uintd<4>::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,4>( typename uintd<4>::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,4>( typename uintd<4>::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,4>( typename uintd<4>::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_expand_with_zero_fill<float,1>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<typename complext<float>::Type,1>( cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_expand_with_zero_fill<float,2>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<typename complext<float>::Type,2>( cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_expand_with_zero_fill<float,3>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<typename complext<float>::Type,3>( cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_expand_with_zero_fill<float,4>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<typename complext<float>::Type,4>( cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,1>(typename uintd<1>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<typename complext<float>::Type,1>(typename uintd<1>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,2>(typename uintd<2>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<typename complext<float>::Type,2>(typename uintd<2>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,3>(typename uintd<3>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<typename complext<float>::Type,3>(typename uintd<3>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,4>(typename uintd<4>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<typename complext<float>::Type,4>(typename uintd<4>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,float,1>(typename reald<float,1>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,typename complext<float>::Type,1>(typename reald<float,1>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,float,2>(typename reald<float,2>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,typename complext<float>::Type,2>(typename reald<float,2>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,float,3>(typename reald<float,3>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,typename complext<float>::Type,3>(typename reald<float,3>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_zero_fill_border<float,float,4>(typename reald<float,4>::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,typename complext<float>::Type,4>(typename reald<float,4>::Type, cuNDArray<typename complext<float>::Type>*);

template void cuNDA_clear<float>(cuNDArray<float>*);
template void cuNDA_clear<typename complext<float>::Type>(cuNDArray<typename complext<float>::Type>*);

template void cuNDA_abs<float>(cuNDArray<float>*);
template void cuNDA_abs<typename floatd<1>::Type>(cuNDArray<typename floatd<1>::Type>*);
template void cuNDA_abs<typename floatd<2>::Type>(cuNDArray<typename floatd<2>::Type>*);
template void cuNDA_abs<typename floatd<3>::Type>(cuNDArray<typename floatd<3>::Type>*);
template void cuNDA_abs<typename floatd<4>::Type>(cuNDArray<typename floatd<4>::Type>*);

template void cuNDA_reciprocal<float>(cuNDArray<float>*);
template void cuNDA_reciprocal<typename complext<float>::Type>(cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_rss_normalize<float, float>(cuNDArray<float>*, unsigned int);
template bool cuNDA_rss_normalize<float, typename complext<float>::Type>(cuNDArray<typename complext<float>::Type>*, unsigned int);

template bool cuNDA_scale<float, float>(cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_scale<float, typename complext<float>::Type>(cuNDArray<float>*, cuNDArray<typename complext<float>::Type>*);
template bool cuNDA_scale<typename complext<float>::Type, typename complext<float>::Type>(cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*);

template void cuNDA_scale<float, float>(float, cuNDArray<float>*);
template void cuNDA_scale<float, typename complext<float>::Type>(float, cuNDArray<typename complext<float>::Type>*);
template void cuNDA_scale<typename complext<float>::Type, typename complext<float>::Type>(typename complext<float>::Type, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_axpy<float, float>( float, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_axpy<float, typename complext<float>::Type>( float, cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*);

template bool cuNDA_axpby<float, float, float>( cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_axpby<float, float, typename complext<float>::Type>( cuNDArray<float>*, cuNDArray<typename complext<float>::Type>*, cuNDArray<float>*, cuNDArray<typename complext<float>::Type>*);
template bool cuNDA_axpby<typename complext<float>::Type, typename complext<float>::Type, typename complext<float>::Type>( cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*, cuNDArray<typename complext<float>::Type>*);

template float cuNDA_dot<float>( cuNDArray<float>*, cuNDArray<float>*, cublasHandle_t );
template float_complext::Type cuNDA_dot<float_complext::Type>( cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*, cublasHandle_t );

template bool cuNDA_axpy<float>( float, cuNDArray<float>*, cuNDArray<float>*, cublasHandle_t );
template bool cuNDA_axpy<typename float_complext::Type>( typename float_complext::Type, cuNDArray<typename float_complext::Type>*, cuNDArray<typename float_complext::Type>*, cublasHandle_t );

template bool cuNDA_scal<float>( float, cuNDArray<float>*, cublasHandle_t );
template bool cuNDA_scal<typename float_complext::Type>( typename float_complext::Type, cuNDArray<typename float_complext::Type>*, cublasHandle_t );

// Instanciation -- double precision

template auto_ptr< cuNDArray<double> > cuNDA_sum<double>( cuNDArray<double>*, unsigned int);
template auto_ptr< cuNDArray<doubled<1>::Type> > cuNDA_sum<doubled<1>::Type>( cuNDArray<doubled<1>::Type>*, unsigned int );
template auto_ptr< cuNDArray<doubled<2>::Type> > cuNDA_sum<doubled<2>::Type>( cuNDArray<doubled<2>::Type>*, unsigned int );
template auto_ptr< cuNDArray<doubled<3>::Type> > cuNDA_sum<doubled<3>::Type>( cuNDArray<doubled<3>::Type>*, unsigned int );
template auto_ptr< cuNDArray<doubled<4>::Type> > cuNDA_sum<doubled<4>::Type>( cuNDArray<doubled<4>::Type>*, unsigned int );

template auto_ptr< cuNDArray<double> > cuNDA_norm<double>( cuNDArray<double>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm<double,1>( cuNDArray<doubled<1>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm<double,2>( cuNDArray<doubled<2>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm<double,3>( cuNDArray<doubled<3>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm<double,4>( cuNDArray<doubled<4>::Type>*);

template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double>( cuNDArray<double>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,1>( cuNDArray<doubled<1>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,2>( cuNDArray<doubled<2>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,3>( cuNDArray<doubled<3>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,4>( cuNDArray<doubled<4>::Type>*);

template auto_ptr< cuNDArray<double> > cuNDA_rss<double, double>( cuNDArray<double>*, unsigned int);
template auto_ptr< cuNDArray<double> > cuNDA_rss<double, typename complext<double>::Type>( cuNDArray<typename complext<double>::Type>*, unsigned int);

template auto_ptr< cuNDArray<typename complext<double>::Type> > cuNDA_correlation<double>( cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_crop<double,1>( typename uintd<1>::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,1>( typename uintd<1>::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,1>( typename uintd<1>::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,1>( typename uintd<1>::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,1>( typename uintd<1>::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_crop<double,2>( typename uintd<2>::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,2>( typename uintd<2>::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,2>( typename uintd<2>::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,2>( typename uintd<2>::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,2>( typename uintd<2>::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_crop<double,3>( typename uintd<3>::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,3>( typename uintd<3>::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,3>( typename uintd<3>::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,3>( typename uintd<3>::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,3>( typename uintd<3>::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_crop<double,4>( typename uintd<4>::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,4>( typename uintd<4>::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,4>( typename uintd<4>::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,4>( typename uintd<4>::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,4>( typename uintd<4>::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_expand_with_zero_fill<double,1>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<typename complext<double>::Type,1>( cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_expand_with_zero_fill<double,2>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<typename complext<double>::Type,2>( cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_expand_with_zero_fill<double,3>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<typename complext<double>::Type,3>( cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_expand_with_zero_fill<double,4>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<typename complext<double>::Type,4>( cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,1>(typename uintd<1>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<typename complext<double>::Type,1>(typename uintd<1>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,2>(typename uintd<2>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<typename complext<double>::Type,2>(typename uintd<2>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,3>(typename uintd<3>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<typename complext<double>::Type,3>(typename uintd<3>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,4>(typename uintd<4>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<typename complext<double>::Type,4>(typename uintd<4>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,double,1>(typename reald<double,1>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,typename complext<double>::Type,1>(typename reald<double,1>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,double,2>(typename reald<double,2>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,typename complext<double>::Type,2>(typename reald<double,2>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,double,3>(typename reald<double,3>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,typename complext<double>::Type,3>(typename reald<double,3>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_zero_fill_border<double,double,4>(typename reald<double,4>::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,typename complext<double>::Type,4>(typename reald<double,4>::Type, cuNDArray<typename complext<double>::Type>*);

template void cuNDA_clear<double>(cuNDArray<double>*);
template void cuNDA_clear<typename complext<double>::Type>(cuNDArray<typename complext<double>::Type>*);

template void cuNDA_abs<double>(cuNDArray<double>*);
template void cuNDA_abs<typename doubled<1>::Type>(cuNDArray<typename doubled<1>::Type>*);
template void cuNDA_abs<typename doubled<2>::Type>(cuNDArray<typename doubled<2>::Type>*);
template void cuNDA_abs<typename doubled<3>::Type>(cuNDArray<typename doubled<3>::Type>*);
template void cuNDA_abs<typename doubled<4>::Type>(cuNDArray<typename doubled<4>::Type>*);

template void cuNDA_reciprocal<double>(cuNDArray<double>*);
template void cuNDA_reciprocal<typename complext<double>::Type>(cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_rss_normalize<double, double>(cuNDArray<double>*, unsigned int);
template bool cuNDA_rss_normalize<double, typename complext<double>::Type>(cuNDArray<typename complext<double>::Type>*, unsigned int);

template void cuNDA_scale<double, double>(double, cuNDArray<double>*);
template void cuNDA_scale<double, typename complext<double>::Type>(double, cuNDArray<typename complext<double>::Type>*);
template void cuNDA_scale<typename complext<double>::Type, typename complext<double>::Type>(typename complext<double>::Type, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_scale<double, double>(cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_scale<double, typename complext<double>::Type>(cuNDArray<double>*, cuNDArray<typename complext<double>::Type>*);
template bool cuNDA_scale<typename complext<double>::Type, typename complext<double>::Type>(cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_axpy<double, double>( double, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_axpy<double, typename complext<double>::Type>( double, cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*);

template bool cuNDA_axpby<double, double, double>( cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_axpby<double, double, typename complext<double>::Type>( cuNDArray<double>*, cuNDArray<typename complext<double>::Type>*, cuNDArray<double>*, cuNDArray<typename complext<double>::Type>*);
template bool cuNDA_axpby<typename complext<double>::Type, typename complext<double>::Type, typename complext<double>::Type>( cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*, cuNDArray<typename complext<double>::Type>*);

template double cuNDA_dot<double>( cuNDArray<double>*, cuNDArray<double>*, cublasHandle_t );
template double_complext::Type cuNDA_dot<double_complext::Type>( cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*, cublasHandle_t );

template bool cuNDA_axpy<double>( double, cuNDArray<double>*, cuNDArray<double>*, cublasHandle_t );
template bool cuNDA_axpy<typename double_complext::Type>( typename double_complext::Type, cuNDArray<typename double_complext::Type>*, cuNDArray<typename double_complext::Type>*, cublasHandle_t );

template bool cuNDA_scal<double>( double, cuNDArray<double>*, cublasHandle_t );
template bool cuNDA_scal<typename double_complext::Type>( typename double_complext::Type, cuNDArray<typename double_complext::Type>*, cublasHandle_t );
