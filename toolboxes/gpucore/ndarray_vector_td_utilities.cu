#include "ndarray_vector_td_utilities.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include "real_utilities.h"
#include "check_CUDA.h"

#include <vector>
#include <cmath>

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
template<class T>  
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
template<class REAL, class T> __global__ 
void cuNDA_norm_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T val = in[idx]; 
    out[idx] = norm<REAL,T>(val);
  }
}

// Norm
template<class REAL, unsigned int D>  
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

// Norm
template<class REAL, class T>  
auto_ptr< cuNDArray<REAL> > cuNDA_norm( cuNDArray<T> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
 
  // Make modulus image
  if( out != 0x0 )
    cuNDA_norm_kernel<REAL,T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<REAL> >(out);
}

// Norm squared
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
template<class REAL, class T> __global__ 
void cuNDA_norm_squared_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T val = in[idx]; 
    out[idx] = norm_squared<REAL,T>(val);
  }
}

// Norm Squared
template<class REAL, unsigned int D>  
auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<typename reald<REAL,D>::Type> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
 
  // Make modulus image
  if( out != 0x0 )
    cuNDA_norm_squared_kernel<REAL,D><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<REAL> >(out);
}

// Norm Squared
template<class REAL, class T>  
auto_ptr< cuNDArray<REAL> > cuNDA_norm_squared( cuNDArray<T> *in )
{
  unsigned int number_of_elements = in->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  cuNDArray<REAL> *out = cuNDArray<REAL>::allocate(in->get_dimensions());
 
  // Make modulus image
  if( out != 0x0 )
    cuNDA_norm_squared_kernel<REAL,T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<REAL> >(out);
}

// SS
template<class REAL, class T> __inline__  __device__ REAL
_ss( unsigned int idx, T *in, unsigned int stride, unsigned int number_of_batches )
{
  unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
  REAL ss = get_zero<REAL>();
  
  for( unsigned int i=0; i<number_of_batches; i++ ) 
    ss += norm_squared<REAL>(in[i*stride+in_idx]);
  
  return ss;
}

// SS
template<class REAL, class T> __global__ void
cuNDA_ss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _ss<REAL,T>(idx, in, stride, number_of_batches); 
  }
}

// SS
template<class REAL, class T>  
auto_ptr< cuNDArray<REAL> > _cuNDA_ss( cuNDArray<T> *in, unsigned int dim )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_ss:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<REAL> >(0x0);
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_ss:: dimension out of range." << endl; 
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
    cuNDA_ss_kernel<REAL,T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<REAL> >(out);
}

template<>  
auto_ptr< cuNDArray<float> > cuNDA_ss<float,float>( cuNDArray<float> *in, unsigned int dim )
{
  return _cuNDA_ss<float, float>(in,dim);
}

template<>  
auto_ptr< cuNDArray<double> > cuNDA_ss<double,double>( cuNDArray<double> *in, unsigned int dim )
{
  return _cuNDA_ss<double, double>(in,dim);
}

template<>
auto_ptr< cuNDArray<float> > cuNDA_ss<float,float_complext::Type>( cuNDArray<float_complext::Type> *in, unsigned int dim )
{
  return _cuNDA_ss<float, float_complext::Type>(in,dim);
}

template<>
auto_ptr< cuNDArray<double> > cuNDA_ss<double,double_complext::Type>( cuNDArray<double_complext::Type> *in, unsigned int dim )
{
  return _cuNDA_ss<double, double_complext::Type>(in,dim);
}

// cSS
template<class REAL, class T> __global__ void
cuNDA_css_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    REAL ss = _ss<REAL,T>(idx, in, stride, number_of_batches); 

    out[idx].vec[0] = ss;
    out[idx].vec[1] = get_zero<REAL>();
  }
}

// cSS
template<class REAL>  
auto_ptr< cuNDArray<typename complext<REAL>::Type> > cuNDA_css( cuNDArray<typename complext<REAL>::Type> *in, unsigned int dim )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_ss:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_ss:: dimension out of range." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
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

  cuNDArray<typename complext<REAL>::Type> *out = cuNDArray<typename complext<REAL>::Type>::allocate(dims);
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  if ( out != 0x0 ){
    cuNDA_css_kernel<REAL, typename complext<REAL>::Type><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  }
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(out);
}

template<>
auto_ptr< cuNDArray<float_complext::Type> > cuNDA_ss<float_complext::Type, float_complext::Type>( cuNDArray<float_complext::Type> *in, unsigned int dim )
{
  return cuNDA_css<float>( in, dim );
}

template<>
auto_ptr< cuNDArray<double_complext::Type> > cuNDA_ss<double_complext::Type, double_complext::Type>( cuNDArray<double_complext::Type> *in, unsigned int dim )
{
  return cuNDA_css<double>( in, dim );
}

// RSS
template<class REAL, class T> __inline__  __device__ REAL
_rss( unsigned int idx, T *in, unsigned int stride, unsigned int number_of_batches )
{
  unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
  REAL rss = get_zero<REAL>();
  
  for( unsigned int i=0; i<number_of_batches; i++ ) 
    rss += norm_squared<REAL>(in[i*stride+in_idx]);
  
  rss = sqrt(rss); 

  return rss;
}

// RSS
template<class REAL, class T> __global__ void
cuNDA_rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _rss<REAL,T>(idx, in, stride, number_of_batches); 
  }
}

// RSS
template<class REAL, class T>  
auto_ptr< cuNDArray<REAL> > _cuNDA_rss( cuNDArray<T> *in, unsigned int dim )
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

template<>  
auto_ptr< cuNDArray<float> > cuNDA_rss<float,float>( cuNDArray<float> *in, unsigned int dim )
{
  return _cuNDA_rss<float, float>(in,dim);
}

template<>  
auto_ptr< cuNDArray<double> > cuNDA_rss<double,double>( cuNDArray<double> *in, unsigned int dim )
{
  return _cuNDA_rss<double, double>(in,dim);
}

template<>
auto_ptr< cuNDArray<float> > cuNDA_rss<float,float_complext::Type>( cuNDArray<float_complext::Type> *in, unsigned int dim )
{
  return _cuNDA_rss<float, float_complext::Type>(in,dim);
}

template<>
auto_ptr< cuNDArray<double> > cuNDA_rss<double,double_complext::Type>( cuNDArray<double_complext::Type> *in, unsigned int dim )
{
  return _cuNDA_rss<double, double_complext::Type>(in,dim);
}

// cRSS
template<class REAL, class T> __global__ void
cuNDA_crss_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    REAL rss = _rss<REAL,T>(idx, in, stride, number_of_batches); 

    out[idx].vec[0] = rss;
    out[idx].vec[1] = get_zero<REAL>();
  }
}

// cRSS
template<class REAL>  
auto_ptr< cuNDArray<typename complext<REAL>::Type> > cuNDA_crss( cuNDArray<typename complext<REAL>::Type> *in, unsigned int dim )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_rss:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_rss:: dimension out of range." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
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

  cuNDArray<typename complext<REAL>::Type> *out = cuNDArray<typename complext<REAL>::Type>::allocate(dims);
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  if ( out != 0x0 ){
    cuNDA_crss_kernel<REAL, typename complext<REAL>::Type><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  }
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(out);
}

template<>
auto_ptr< cuNDArray<float_complext::Type> > cuNDA_rss<float_complext::Type, float_complext::Type>( cuNDArray<float_complext::Type> *in, unsigned int dim )
{
  return cuNDA_crss<float>( in, dim );
}

template<>
auto_ptr< cuNDArray<double_complext::Type> > cuNDA_rss<double_complext::Type, double_complext::Type>( cuNDArray<double_complext::Type> *in, unsigned int dim )
{
  return cuNDA_crss<double>( in, dim );
}

// RSS
template<class REAL, class T> __global__ void
cuNDA_reciprocal_rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = reciprocal<REAL>(_rss<REAL,T>(idx, in, stride, number_of_batches));
  }
}

// Reciprocal RSS
template<class REAL, class T>  
auto_ptr< cuNDArray<REAL> > _cuNDA_reciprocal_rss( cuNDArray<T> *in, unsigned int dim )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_reciprocal_rss:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<REAL> >(0x0);
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_reciprocal_rss:: dimension out of range." << endl; 
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
    cuNDA_reciprocal_rss_kernel<REAL,T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<REAL> >(out);
}

template<>
auto_ptr< cuNDArray<float> > cuNDA_reciprocal_rss<float,float>( cuNDArray<float> *in, unsigned int dim )
{
  return _cuNDA_reciprocal_rss<float, float>(in,dim);
}

template<>
auto_ptr< cuNDArray<double> > cuNDA_reciprocal_rss<double,double>( cuNDArray<double> *in, unsigned int dim )
{
  return _cuNDA_reciprocal_rss<double, double>(in,dim);
}

template<>
auto_ptr< cuNDArray<float> > cuNDA_reciprocal_rss<float,float_complext::Type>( cuNDArray<float_complext::Type> *in, unsigned int dim )
{
  return _cuNDA_reciprocal_rss<float, float_complext::Type>(in,dim);
}

template<>
auto_ptr< cuNDArray<double> > cuNDA_reciprocal_rss<double,double_complext::Type>( cuNDArray<double_complext::Type> *in, unsigned int dim )
{
  return _cuNDA_reciprocal_rss<double, double_complext::Type>(in,dim);
}

// cReciprocal RSS
template<class REAL, class T> __global__ void
cuNDA_creciprocal_rss_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    REAL reciprocal_rss = reciprocal<REAL>(_rss<REAL,T>(idx, in, stride, number_of_batches));

    out[idx].vec[0] = reciprocal_rss;
    out[idx].vec[1] = get_zero<REAL>();
  }
}

// cReciprocal RSS
template<class REAL>  
auto_ptr< cuNDArray<typename complext<REAL>::Type> > cuNDA_creciprocal_rss( cuNDArray<typename complext<REAL>::Type> *in, unsigned int dim )
{
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_reciprocal_rss:: underdimensioned." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_reciprocal_rss:: dimension out of range." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
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

  cuNDArray<typename complext<REAL>::Type> *out = cuNDArray<typename complext<REAL>::Type>::allocate(dims);
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if ( out != 0x0 )
    cuNDA_creciprocal_rss_kernel<REAL, typename complext<REAL>::Type><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(out);
}

template<>
auto_ptr< cuNDArray<float_complext::Type> > cuNDA_reciprocal_rss<float_complext::Type, float_complext::Type>( cuNDArray<float_complext::Type> *in, unsigned int dim )
{
  return cuNDA_creciprocal_rss<float>( in, dim );
}

template<>
auto_ptr< cuNDArray<double_complext::Type> > cuNDA_reciprocal_rss<double_complext::Type, double_complext::Type>( cuNDArray<double_complext::Type> *in, unsigned int dim )
{
  return cuNDA_creciprocal_rss<double>( in, dim );
}

// Build correlation matrix
template<class REAL, class T> __global__ void
cuNDA_correlation_kernel( T *in, T *corrm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.y;

  if( p < num_elements ){
    for( unsigned int j=0; j<i; j++){
      T tmp = mul<T>(in[i*num_elements+p], conj<T>(in[j*num_elements+p]));
      corrm[(j*num_batches+i)*num_elements+p] = tmp;
      corrm[(i*num_batches+j)*num_elements+p] = conj<T>(tmp);
    }
    T tmp = in[i*num_elements+p];
    corrm[(i*num_batches+i)*num_elements+p] = mul<T>(tmp,conj<T>(tmp));
  }
}

// Build correlation matrix
template<class REAL, class T>  
auto_ptr< cuNDArray<T> > _cuNDA_correlation( cuNDArray<T> *in )
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

  dim3 blockDim(((deviceProp.maxThreadsPerBlock/number_of_batches)/warp_size)*warp_size, number_of_batches);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  if( blockDim.x == 0 ){
    cout << endl << "cuNDA_correlation:: correlation dimension exceeds capacity." << endl; 
    return auto_ptr< cuNDArray<T> >(0x0);
  }

  if( out != 0x0 )
    cuNDA_correlation_kernel<REAL,T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  return auto_ptr< cuNDArray<T> >(out);
}

// Build correlation matrix
template<class REAL> __global__ void
cuNDA_real_to_complext_kernel( REAL *in, typename complext<REAL>::Type *out, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){
    typename complext<REAL>::Type z;
    z.vec[0] = in[idx];
    z.vec[1] = get_zero<REAL>();
    out[idx] = z;
  }
}

// Build correlation matrix
template<class REAL>  
auto_ptr< cuNDArray<typename complext<REAL>::Type> > cuNDA_real_to_complext( cuNDArray<REAL> *in )
{
  cuNDArray<typename complext<REAL>::Type> *out = cuNDArray<typename complext<REAL>::Type>::allocate(in->get_dimensions());

  dim3 blockDim( 256 );
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x));
  
  if( out != 0x0 )
    cuNDA_real_to_complext_kernel<REAL><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements());
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(out);
}

template<>
auto_ptr< cuNDArray<float> > cuNDA_correlation<float>( cuNDArray<float> *data )
{
  return _cuNDA_correlation<float,float>(data);
}

template<>
auto_ptr< cuNDArray<double> > cuNDA_correlation<double>( cuNDArray<double> *data )
{
  return _cuNDA_correlation<double,double>(data);
}

template<>
auto_ptr< cuNDArray<float_complext::Type> > cuNDA_correlation<float_complext::Type>( cuNDArray<float_complext::Type> *data )
{
  return _cuNDA_correlation<float,float_complext::Type>(data);
}

template<>
auto_ptr< cuNDArray<double_complext::Type> > cuNDA_correlation<double_complext::Type>( cuNDArray<double_complext::Type> *data )
{
  return _cuNDA_correlation<double,double_complext::Type>(data);
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
template<class T> 
void cuNDA_clear( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make clear image
  cuNDA_clear_kernel<T><<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
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
template<class T>  
void cuNDA_abs( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
 
  // Make modulus image
  cuNDA_abs_kernel<T><<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
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
template<class T> 
void cuNDA_reciprocal( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_reciprocal_kernel<T><<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
}

// Square root
template<class T> __global__ 
void cuNDA_sqrt_kernel( T *in_out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = sqrt(in_out[idx]);
  }
}

// Square root
template<class T> 
void cuNDA_sqrt( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_sqrt_kernel<T><<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
}

// Reciprocal square root
template<class T> __global__ 
void cuNDA_reciprocal_sqrt_kernel( T *in_out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = gad_rsqrt(in_out[idx]);
  }
}

// Square root
template<class T> 
void cuNDA_reciprocal_sqrt( cuNDArray<T> *in_out )
{
  unsigned int number_of_elements = in_out->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Make reciprocal image
  cuNDA_reciprocal_sqrt_kernel<T><<< gridDim, blockDim >>>( in_out->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
}

// Normalize (float)
template<> void cuNDA_normalize<float>( cuNDArray<float> *data, float new_max, cublasHandle_t handle )
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
template<> void cuNDA_normalize<double>( cuNDArray<double> *data, double new_max, cublasHandle_t handle )
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

// Normalized RSS
template<class REAL, class T> __global__ void
cuNDA_rss_normalize_kernel( T *in_out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < number_of_elements ){

    REAL reciprocal_rss = reciprocal<REAL>(_rss<REAL,T>(idx, in_out, stride, number_of_batches));
 
    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

    for( unsigned int i=0; i<number_of_batches; i++ ) {
      T out = in_out[i*stride+in_idx];
      out *= reciprocal_rss; // complex-scalar multiplication (element-wise operator)
      in_out[i*stride+in_idx] = out; 
    } 
  }
}

// Normalized RSS
template<class REAL, class T> 
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
template<class REAL> __global__ 
void cuNDA_scale1_kernel( REAL a, typename complext<REAL>::Type *x, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    typename complext<REAL>::Type in = x[idx];
    in = mul<REAL>(a,in);
    x[idx] = in;
  }
}

// Scale 
template<class REAL> 
void cuNDA_scale( REAL a, cuNDArray<typename complext<REAL>::Type> *x )
{
  unsigned int number_of_elements = x->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_scale1_kernel<REAL><<< gridDim, blockDim >>> ( a, x->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
}

// Scale
template<class S, class T> __global__ 
void cuNDA_scale2_kernel( S *a, T *x, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    unsigned int frame_offset = blockIdx.y*number_of_elements;
    S in_a = a[idx];
    T in_x = x[idx+frame_offset];
    x[idx+frame_offset] = mul<S>(in_a,in_x);
  }
}

// Scale 
template<class T> 
bool cuNDA_scale( cuNDArray<T> *a, cuNDArray<T> *x )
{
  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch, cannot scale" << endl;
    return false;
  }
 
  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x), num_batches);
 
  // Invoke kernel
  cuNDA_scale2_kernel<T,T><<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// Scale 
template<class REAL> 
bool cuNDA_scale( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x )
{
  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch, cannot scale" << endl;
    return false;
  }
 
  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x), num_batches);
 
  // Invoke kernel
  cuNDA_scale2_kernel<REAL, typename complext<REAL>::Type><<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// 'axpy'
template<class S, class T> __global__ 
void cuNDA_axpy_kernel( S *a, T *x, T *y, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < number_of_elements ){
    S in_a = a[idx];
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      unsigned int iidx = batch*number_of_elements + idx;
      T in_x = x[iidx];
      T in_y = y[iidx];
      in_y += mul<S>(in_a,in_x);
      y[iidx] = in_y;
    }
  }
}

// '.axpy' 
template<class T> 
bool cuNDA_axpy( cuNDArray<T> *a, cuNDArray<T> *x, cuNDArray<T> *y )
{
  if( x->get_number_of_elements() != y->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpy'" << endl;
    return false;
  }

  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpy'" << endl;
    return false;
  }
 
  unsigned int number_of_batches = x->get_number_of_elements() / a->get_number_of_elements();
  unsigned int number_of_elements = a->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_axpy_kernel<T,T><<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), y->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// '.axpy' 
template<class REAL> 
bool cuNDA_axpy( cuNDArray<REAL> *a, cuNDArray<typename complext<REAL>::Type> *x, cuNDArray<typename complext<REAL>::Type> *y )
{
  if( x->get_number_of_elements() != y->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpy'" << endl;
    return false;
  }

  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch in 'axpy'" << endl;
    return false;
  }
 
  unsigned int number_of_batches = x->get_number_of_elements() / a->get_number_of_elements();
  unsigned int number_of_elements = a->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  // Invoke kernel
  cuNDA_axpy_kernel<REAL, typename complext<REAL>::Type><<< gridDim, blockDim >>> ( a->get_data_ptr(), x->get_data_ptr(), y->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
  return true;
}

// Crop
template<class T, unsigned int D> __global__ void
cuNDA_crop_kernel( typename uintd<D>::Type offset, typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, T *in, T *out )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < prod(matrix_size_out) ){
    const typename uintd<D>::Type co = idx_to_co<D>( idx, matrix_size_out );
    const typename uintd<D>::Type co_os = offset + co;
    const unsigned int in_idx = co_to_idx<D>(co_os, matrix_size_in)+blockIdx.y*prod(matrix_size_in);
    const unsigned int out_idx = idx+blockIdx.y*prod(matrix_size_out);
    out[out_idx] = in[in_idx];
  }
}

// Crop
template<class T, unsigned int D> 
bool cuNDA_crop( typename uintd<D>::Type offset, cuNDArray<T> *in, cuNDArray<T> *out )
{ 
  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "image dimensions mismatch, cannot crop" << endl;
    return false;
  }

  if( in->get_number_of_dimensions() < D ){
    cout << endl << "number of image dimensions should be at least " << D << ", cannot crop" << endl;
    return false;
  }

  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( out->get_dimensions() );
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ )
    number_of_batches *= in->get_size(d);

  if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
    cout << endl << "cropping size mismatch, cannot crop" << endl;
    return false;
  }

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)prod(matrix_size_out)/blockDim.x), number_of_batches);

  // Invoke kernel
  cuNDA_crop_kernel<T,D><<< gridDim, blockDim >>> ( offset, matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr() );
 
  CHECK_FOR_CUDA_ERROR();

  return true;
}

// Expand and zero fill
template<class T, unsigned int D> __global__ void
cuNDA_expand_with_zero_fill_kernel( typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, T *in, T *out )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < prod(matrix_size_out) ){

    const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
    const typename uintd<D>::Type offset = (matrix_size_out-matrix_size_in)>>1;
    T _out;
    bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));

    if( inside )
      _out = in[co_to_idx<D>(co_out-offset, matrix_size_in)+blockIdx.y*prod(matrix_size_in)];
    else{      
      _out = get_zero<T>();
    }

    out[idx+blockIdx.y*prod(matrix_size_out)] = _out;
  }
}

// Expand and zero fill
template<class T, unsigned int D> 
bool cuNDA_expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out )
{ 
  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "Image dimensions mismatch, cannot expand" << endl;
    return false;
  }
 
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( out->get_dimensions() );
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ )
    number_of_batches *= in->get_size(d);

  if( weak_greater(matrix_size_in,matrix_size_out) ){
    cout << endl << "Size mismatch, cannot expand" << endl;
    return false;
  }
 
  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)prod(matrix_size_out)/blockDim.x), number_of_batches );
 
  // Invoke kernel
  cuNDA_expand_with_zero_fill_kernel<T,D><<< gridDim, blockDim >>> 
    ( matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr() );
 
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
template<class T, unsigned int D> 
bool cuNDA_zero_fill_border( typename uintd<D>::Type matrix_size_in, cuNDArray<T> *out )
{ 
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( out->get_dimensions() );
 
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
    
    typename reald<REAL,D>::Type half_matrix_size_out_real = to_reald<REAL,unsigned int,D>( matrix_size_out>>1 );

    const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
    typename reald<REAL,D>::Type co_out_real = to_reald<REAL,unsigned int,D>( co_out );
    
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
template<class REAL, class T, unsigned int D> 
bool cuNDA_zero_fill_border( typename reald<REAL,D>::Type radius, cuNDArray<T> *out )
{
  if( out->get_number_of_dimensions() != D ){
    cout << endl << "Image dimensions mismatch, cannot zero fill" << endl;
    return false;
  }
 
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( out->get_dimensions() );
  typename reald<REAL,D>::Type matrix_size_out_real = to_reald<REAL,unsigned int,D>( matrix_size_out );

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

template<> float_complext::Type
_dot<float_complext::Type>( cuNDArray<float_complext::Type>* arr1, cuNDArray<float_complext::Type>* arr2, cublasHandle_t handle )
{
  float_complext::Type ret;
  
  if (cublasCdotc( handle, arr1->get_number_of_elements(),
		   (const cuComplex*) arr1->get_data_ptr(), 1, 
		   (const cuComplex*) arr2->get_data_ptr(), 1,
		   (cuComplex*) &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> float 
_dot<float>( cuNDArray<float>* arr1, cuNDArray<float>* arr2, cublasHandle_t handle )
{
  float ret;
  if( cublasSdot(handle, arr1->get_number_of_elements(),
		 arr1->get_data_ptr(), 1, 
		 arr2->get_data_ptr(), 1,
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> double_complext::Type
_dot<double_complext::Type>( cuNDArray<double_complext::Type>* arr1, cuNDArray<double_complext::Type>* arr2, cublasHandle_t handle )
{
  double_complext::Type ret;
  
  if( cublasZdotc(handle, arr1->get_number_of_elements(),
		  (const cuDoubleComplex*) arr1->get_data_ptr(), 1, 
		  (const cuDoubleComplex*) arr2->get_data_ptr(), 1,
		  (cuDoubleComplex*) &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> double 
_dot<double>( cuNDArray<double>* arr1, cuNDArray<double>* arr2, cublasHandle_t handle )
{
  double ret;
  if( cublasDdot(handle, arr1->get_number_of_elements(),
		 arr1->get_data_ptr(), 1, 
		 arr2->get_data_ptr(), 1,
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
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

template<class REAL, class T> REAL _sum( cuNDArray<T>* arr, cublasHandle_t handle );

template<> float
_sum<float,float_complext::Type>( cuNDArray<float_complext::Type>* arr, cublasHandle_t handle )
{
  float ret;
  
  if (cublasScasum( handle, arr->get_number_of_elements(),
		   (const cuComplex*) arr->get_data_ptr(), 1,
		   &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> float 
_sum<float,float>( cuNDArray<float>* arr, cublasHandle_t handle )
{
  float ret;
  if( cublasSasum(handle, arr->get_number_of_elements(),
		 arr->get_data_ptr(), 1, 
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> double
_sum<double,double_complext::Type>( cuNDArray<double_complext::Type>* arr, cublasHandle_t handle )
{
  double ret;
  
  if( cublasDzasum(handle, arr->get_number_of_elements(),
		  (const cuDoubleComplex*) arr->get_data_ptr(), 1, 
		  &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
    }
  
  return ret;
}

template<> double 
_sum<double,double>( cuNDArray<double>* arr, cublasHandle_t handle )
{
  double ret;
  if( cublasDasum(handle, arr->get_number_of_elements(),
		 arr->get_data_ptr(), 1, 
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
    }
  
  return ret;
}

template<class REAL, class T> REAL
cuNDA_asum( cuNDArray<T>* arr, cublasHandle_t handle )
{
  return _sum<REAL,T>( arr, handle );  
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
      cout << "cuNDA_axpy: axpy calculating using cublas failed" << std::endl;
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
      cout << "cuNDA_axpy: axpy calculating using cublas failed" << std::endl;
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
      cout << "cuNDA_axpy: axpy calculating using cublas failed" << std::endl;
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
      cout << "cuNDA_axpy: axpy calculating using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<class T> bool
cuNDA_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y, cublasHandle_t handle )
{
  if (x->get_number_of_elements() != y->get_number_of_elements()) {
    cout << "cuNDA_axpy: axpy array dimensions mismatch" << std::endl;
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

//
// Instantiation
//

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

template auto_ptr< cuNDArray<float> > cuNDA_norm<float,float>( cuNDArray<float>*);

template auto_ptr< cuNDArray<float> > cuNDA_norm<float,float_complext::Type>( cuNDArray<float_complext::Type>* );

template auto_ptr< cuNDArray<float> > cuNDA_norm<float,1>( cuNDArray<floatd<1>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm<float,2>( cuNDArray<floatd<2>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm<float,3>( cuNDArray<floatd<3>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm<float,4>( cuNDArray<floatd<4>::Type>*);

template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,float>( cuNDArray<float>*);

template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,float_complext::Type>( cuNDArray<float_complext::Type>* );

template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,1>( cuNDArray<floatd<1>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,2>( cuNDArray<floatd<2>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,3>( cuNDArray<floatd<3>::Type>*);
template auto_ptr< cuNDArray<float> > cuNDA_norm_squared<float,4>( cuNDArray<floatd<4>::Type>*);

template auto_ptr< cuNDArray<float> > cuNDA_ss<float,float>( cuNDArray<float>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_ss<float,float_complext::Type>( cuNDArray<float_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<float_complext::Type> > cuNDA_ss<float_complext::Type, float_complext::Type>( cuNDArray<float_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_rss<float,float>( cuNDArray<float>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_rss<float,float_complext::Type>( cuNDArray<float_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<float_complext::Type> > cuNDA_rss<float_complext::Type, float_complext::Type>( cuNDArray<float_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_reciprocal_rss<float,float>( cuNDArray<float>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_reciprocal_rss<float,float_complext::Type>( cuNDArray<float_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<float_complext::Type> > cuNDA_reciprocal_rss<float_complext::Type, float_complext::Type>( cuNDArray<float_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<float> > cuNDA_correlation<float>( cuNDArray<float>*);
template auto_ptr< cuNDArray<float_complext::Type> > cuNDA_correlation<float_complext::Type>( cuNDArray<float_complext::Type>*);

template auto_ptr< cuNDArray<float_complext::Type> > cuNDA_real_to_complext<float>( cuNDArray<float>*);

template void cuNDA_clear<float>(cuNDArray<float>*);
template void cuNDA_clear<float_complext::Type>(cuNDArray<float_complext::Type>*);

template void cuNDA_reciprocal<float>(cuNDArray<float>*);
template void cuNDA_reciprocal<float_complext::Type>(cuNDArray<float_complext::Type>*);

template void cuNDA_sqrt<float>(cuNDArray<float>*);

template void cuNDA_reciprocal_sqrt<float>(cuNDArray<float>*);

template void cuNDA_abs<float>(cuNDArray<float>*);
template void cuNDA_abs<floatd1::Type>(cuNDArray<floatd1::Type>*);
template void cuNDA_abs<floatd2::Type>(cuNDArray<floatd2::Type>*);
template void cuNDA_abs<floatd3::Type>(cuNDArray<floatd3::Type>*);
template void cuNDA_abs<floatd4::Type>(cuNDArray<floatd4::Type>*);

template bool cuNDA_rss_normalize<float,float>(cuNDArray<float>*, unsigned int);
template bool cuNDA_rss_normalize<float,float_complext::Type>(cuNDArray<float_complext::Type>*, unsigned int);

template void cuNDA_scale<float>(float, cuNDArray<float_complext::Type>*);

template bool cuNDA_scale<float>(cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_scale<float_complext::Type>(cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_scale<float>(cuNDArray<float>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_axpy<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_axpy<float_complext::Type>( cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_axpy<float>( cuNDArray<float>*, cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_crop<float,1>( uintd1::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,1>( uintd1::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,1>( uintd1::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,1>( uintd1::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,1>( uintd1::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_crop<float,2>( uintd2::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,2>( uintd2::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,2>( uintd2::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,2>( uintd2::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,2>( uintd2::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_crop<float,3>( uintd3::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,3>( uintd3::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,3>( uintd3::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,3>( uintd3::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,3>( uintd3::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_crop<float,4>( uintd4::Type, cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_crop<vector_td<float,1>,4>( uintd4::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*);
template bool cuNDA_crop<vector_td<float,2>,4>( uintd4::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*);
template bool cuNDA_crop<vector_td<float,3>,4>( uintd4::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*);
template bool cuNDA_crop<vector_td<float,4>,4>( uintd4::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*);

template bool cuNDA_expand_with_zero_fill<float,1>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<float_complext::Type,1>( cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_expand_with_zero_fill<float,2>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<float_complext::Type,2>( cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_expand_with_zero_fill<float,3>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<float_complext::Type,3>( cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_expand_with_zero_fill<float,4>( cuNDArray<float>*, cuNDArray<float>*);
template bool cuNDA_expand_with_zero_fill<float_complext::Type,4>( cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,1>(uintd1::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float_complext::Type,1>(uintd1::Type, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,2>(uintd2::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float_complext::Type,2>(uintd2::Type, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,3>(uintd3::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float_complext::Type,3>(uintd3::Type, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,4>(uintd4::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float_complext::Type,4>(uintd4::Type, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,float,1>(floatd1::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,float_complext::Type,1>(floatd1::Type, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,float,2>(floatd2::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,float_complext::Type,2>(floatd2::Type, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,float,3>(floatd3::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,float_complext::Type,3>(floatd3::Type, cuNDArray<float_complext::Type>*);

template bool cuNDA_zero_fill_border<float,float,4>(floatd4::Type, cuNDArray<float>*);
template bool cuNDA_zero_fill_border<float,float_complext::Type,4>(floatd4::Type, cuNDArray<float_complext::Type>*);


template float cuNDA_dot<float>( cuNDArray<float>*, cuNDArray<float>*, cublasHandle_t );
template float_complext::Type cuNDA_dot<float_complext::Type>( cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*, cublasHandle_t );

template float cuNDA_asum<float,float>( cuNDArray<float>*, cublasHandle_t );
template float cuNDA_asum<float,float_complext::Type>( cuNDArray<float_complext::Type>*, cublasHandle_t );

template bool cuNDA_axpy<float>( float, cuNDArray<float>*, cuNDArray<float>*, cublasHandle_t );
template bool cuNDA_axpy<float_complext::Type>( float_complext::Type, cuNDArray<float_complext::Type>*, cuNDArray<float_complext::Type>*, cublasHandle_t );

template bool cuNDA_scal<float>( float, cuNDArray<float>*, cublasHandle_t );
template bool cuNDA_scal<float_complext::Type>( float_complext::Type, cuNDArray<float_complext::Type>*, cublasHandle_t );

// Instanciation -- double precision

template auto_ptr< cuNDArray<double> > cuNDA_sum<double>( cuNDArray<double>*, unsigned int);
template auto_ptr< cuNDArray<doubled<1>::Type> > cuNDA_sum<doubled<1>::Type>( cuNDArray<doubled<1>::Type>*, unsigned int );
template auto_ptr< cuNDArray<doubled<2>::Type> > cuNDA_sum<doubled<2>::Type>( cuNDArray<doubled<2>::Type>*, unsigned int );
template auto_ptr< cuNDArray<doubled<3>::Type> > cuNDA_sum<doubled<3>::Type>( cuNDArray<doubled<3>::Type>*, unsigned int );
template auto_ptr< cuNDArray<doubled<4>::Type> > cuNDA_sum<doubled<4>::Type>( cuNDArray<doubled<4>::Type>*, unsigned int );

template auto_ptr< cuNDArray<double> > cuNDA_norm<double,double>( cuNDArray<double>*);

template auto_ptr< cuNDArray<double> > cuNDA_norm<double,double_complext::Type>( cuNDArray<double_complext::Type>* );

template auto_ptr< cuNDArray<double> > cuNDA_norm<double,1>( cuNDArray<doubled<1>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm<double,2>( cuNDArray<doubled<2>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm<double,3>( cuNDArray<doubled<3>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm<double,4>( cuNDArray<doubled<4>::Type>*);

template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,double>( cuNDArray<double>*);

template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,double_complext::Type>( cuNDArray<double_complext::Type>* );

template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,1>( cuNDArray<doubled<1>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,2>( cuNDArray<doubled<2>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,3>( cuNDArray<doubled<3>::Type>*);
template auto_ptr< cuNDArray<double> > cuNDA_norm_squared<double,4>( cuNDArray<doubled<4>::Type>*);

template auto_ptr< cuNDArray<double> > cuNDA_ss<double,double>( cuNDArray<double>*, unsigned int);

template auto_ptr< cuNDArray<double> > cuNDA_ss<double,double_complext::Type>( cuNDArray<double_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<double_complext::Type> > cuNDA_ss<double_complext::Type, double_complext::Type>( cuNDArray<double_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<double> > cuNDA_rss<double,double>( cuNDArray<double>*, unsigned int);

template auto_ptr< cuNDArray<double> > cuNDA_rss<double,double_complext::Type>( cuNDArray<double_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<double_complext::Type> > cuNDA_rss<double_complext::Type, double_complext::Type>( cuNDArray<double_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<double> > cuNDA_reciprocal_rss<double,double>( cuNDArray<double>*, unsigned int);

template auto_ptr< cuNDArray<double> > cuNDA_reciprocal_rss<double,double_complext::Type>( cuNDArray<double_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<double_complext::Type> > cuNDA_reciprocal_rss<double_complext::Type, double_complext::Type>( cuNDArray<double_complext::Type>*, unsigned int);

template auto_ptr< cuNDArray<double> > cuNDA_correlation<double>( cuNDArray<double>*);
template auto_ptr< cuNDArray<double_complext::Type> > cuNDA_correlation<double_complext::Type>( cuNDArray<double_complext::Type>*);

template auto_ptr< cuNDArray<double_complext::Type> > cuNDA_real_to_complext<double>( cuNDArray<double>*);

template void cuNDA_clear<double>(cuNDArray<double>*);
template void cuNDA_clear<double_complext::Type>(cuNDArray<double_complext::Type>*);

template void cuNDA_reciprocal<double>(cuNDArray<double>*);
template void cuNDA_reciprocal<double_complext::Type>(cuNDArray<double_complext::Type>*);

template void cuNDA_sqrt<double>(cuNDArray<double>*);

template void cuNDA_reciprocal_sqrt<double>(cuNDArray<double>*);

template void cuNDA_abs<double>(cuNDArray<double>*);
template void cuNDA_abs<doubled1::Type>(cuNDArray<doubled1::Type>*);
template void cuNDA_abs<doubled2::Type>(cuNDArray<doubled2::Type>*);
template void cuNDA_abs<doubled3::Type>(cuNDArray<doubled3::Type>*);
template void cuNDA_abs<doubled4::Type>(cuNDArray<doubled4::Type>*);

template bool cuNDA_rss_normalize<double,double>(cuNDArray<double>*, unsigned int);
template bool cuNDA_rss_normalize<double,double_complext::Type>(cuNDArray<double_complext::Type>*, unsigned int);

template void cuNDA_scale<double>(double, cuNDArray<double_complext::Type>*);

template bool cuNDA_scale<double>(cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_scale<double_complext::Type>(cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_scale<double>(cuNDArray<double>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_axpy<double>( cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_axpy<double_complext::Type>( cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_axpy<double>( cuNDArray<double>*, cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_crop<double,1>( uintd1::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,1>( uintd1::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,1>( uintd1::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,1>( uintd1::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,1>( uintd1::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_crop<double,2>( uintd2::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,2>( uintd2::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,2>( uintd2::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,2>( uintd2::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,2>( uintd2::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_crop<double,3>( uintd3::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,3>( uintd3::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,3>( uintd3::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,3>( uintd3::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,3>( uintd3::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_crop<double,4>( uintd4::Type, cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_crop<vector_td<double,1>,4>( uintd4::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*);
template bool cuNDA_crop<vector_td<double,2>,4>( uintd4::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*);
template bool cuNDA_crop<vector_td<double,3>,4>( uintd4::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*);
template bool cuNDA_crop<vector_td<double,4>,4>( uintd4::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*);

template bool cuNDA_expand_with_zero_fill<double,1>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<double_complext::Type,1>( cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_expand_with_zero_fill<double,2>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<double_complext::Type,2>( cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_expand_with_zero_fill<double,3>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<double_complext::Type,3>( cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_expand_with_zero_fill<double,4>( cuNDArray<double>*, cuNDArray<double>*);
template bool cuNDA_expand_with_zero_fill<double_complext::Type,4>( cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,1>(uintd1::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double_complext::Type,1>(uintd1::Type, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,2>(uintd2::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double_complext::Type,2>(uintd2::Type, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,3>(uintd3::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double_complext::Type,3>(uintd3::Type, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,4>(uintd4::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double_complext::Type,4>(uintd4::Type, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,double,1>(doubled1::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,double_complext::Type,1>(doubled1::Type, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,double,2>(doubled2::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,double_complext::Type,2>(doubled2::Type, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,double,3>(doubled3::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,double_complext::Type,3>(doubled3::Type, cuNDArray<double_complext::Type>*);

template bool cuNDA_zero_fill_border<double,double,4>(doubled4::Type, cuNDArray<double>*);
template bool cuNDA_zero_fill_border<double,double_complext::Type,4>(doubled4::Type, cuNDArray<double_complext::Type>*);


template double cuNDA_dot<double>( cuNDArray<double>*, cuNDArray<double>*, cublasHandle_t );
template double_complext::Type cuNDA_dot<double_complext::Type>( cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*, cublasHandle_t );

template double cuNDA_asum<double,double>( cuNDArray<double>*, cublasHandle_t );
template double cuNDA_asum<double,double_complext::Type>( cuNDArray<double_complext::Type>*, cublasHandle_t );

template bool cuNDA_axpy<double>( double, cuNDArray<double>*, cuNDArray<double>*, cublasHandle_t );
template bool cuNDA_axpy<double_complext::Type>( double_complext::Type, cuNDArray<double_complext::Type>*, cuNDArray<double_complext::Type>*, cublasHandle_t );

template bool cuNDA_scal<double>( double, cuNDArray<double>*, cublasHandle_t );
template bool cuNDA_scal<double_complext::Type>( double_complext::Type, cuNDArray<double_complext::Type>*, cublasHandle_t );
