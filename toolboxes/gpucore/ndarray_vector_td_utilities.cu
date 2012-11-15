#include "ndarray_vector_td_utilities.h"
#include "real_utilities.h"
#include "real_utilities_device.h"
#include "check_CUDA.h"

#include <cublas_v2.h>

#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

// Some device properties we query once to eliminate runtime overhead
//
static int num_devices = 0;
static int *warp_size = 0x0;
static int *max_blockdim = 0x0;
static int *max_griddim = 0x0;
static cublasHandle_t *handle = 0x0;

// Default template arguments seems to require c++-0x, which we can't assume. 
// We use a dummy type instead...
typedef float dummy;

//
// Some internal utilities
//

// Prepare array for processing:
//
// Sets the device context to the denoted compute device and makes device copies if necessary.
// If compute mode is CUNDA_NDARRAY_DEVICE the mandatory array 'in1' determines the compute device
//

template< unsigned int D, typename I1, typename I2, typename I3 > 
static void prepare( int compute_device, int *cur_device, int *old_device,
		     cuNDArray<I1> *in1,       cuNDArray<I1> **in1_int, 
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> **in2_int = 0x0, 
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> **in3_int = 0x0 ) 
{
  // Test validity of D
  if( D==0 || D>3 ){
    throw std::runtime_error( ">>>Internal error<<< :prepare: D out of range");

  }

  if( !cur_device || !old_device ){
    throw std::runtime_error( ">>>Internal error<<< :prepare: device ids 0x0");

  }

  // Test validity of input pointer
  if( !in1 || !in1_int ){
    throw std::runtime_error( "unable to process 0x0 input");

  }
  if( D>1 && (!in2 || !in2_int) ){
    throw std::runtime_error( "unable to process 0x0 input");

  }
  if( D>2 && (!in3 || !in3_int) ){
    throw std::runtime_error( "unable to process 0x0 input");

  }
  
  // Get current Cuda device
  if( cudaGetDevice(old_device) != cudaSuccess ) {
    throw std::runtime_error( "unable to get device no");

  }

  // Set the cuda device to use for computation
  if( compute_device == CUNDA_CURRENT_DEVICE ){
    *cur_device = *old_device; 
  }
  else if( compute_device == CUNDA_NDARRAY_DEVICE ){
    // Let D indicate which ndarray that determines the device
    // D denotes the output array (the latter ndarray in the list), if any, otherwise a sole input ndarray
    if( D == 1)
      *cur_device = in1->get_device();
    else if( D == 2 )
      *cur_device = in2->get_device();
    else if( D == 3 )
      *cur_device = in3->get_device();
  }
  else{
    throw std::runtime_error( ">>>Internal error<<< :prepare: unknown compute mode");

  }

  if( *cur_device != *old_device && cudaSetDevice(*cur_device) != cudaSuccess) {
    throw std::runtime_error( "unable to set device no");

  }
  
  // Transfer arrays to compute device if necessary
  if( *cur_device != in1->get_device() )
    *in1_int = new cuNDArray<I1>(*in1); // device transfer
  else
    *in1_int = in1;
  
  if( D>1 ){
    if( *cur_device != in2->get_device() )
      *in2_int = new cuNDArray<I2>(*in2); // device transfer
    else
      *in2_int = in2;
  }
  
  if( D>2 ){
    if( *cur_device != in3->get_device() )
      *in3_int = new cuNDArray<I3>(*in3); // device transfer
    else
      *in3_int = in3;
  }
  

}  

// Restore active device and free internal memory 
//
template< unsigned int D, typename I1, typename O, typename I2, typename I3 > 
static void restore( int old_device,
		     cuNDArray<I1> *in1, cuNDArray<I1> *in1_int, 
		     unsigned int out_idx = 0, cuNDA_device alloc_device = CUNDA_NDARRAY_DEVICE, cuNDArray<O>  *out = 0x0,
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> *in2_int = 0x0, 
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> *in3_int = 0x0 )
{
  // Test validity of D
  if( D==0 || D>3 ){
    throw std::runtime_error( ">>>Internal error<<< :prepare: D out of range");

  }

  // Test validity of input pointer
  if( !in1 || !in1_int ){
    throw std::runtime_error( "unable to process 0x0 input");

  }
  if( D>1 && (!in2 || !in2_int) ){
    throw std::runtime_error( "unable to process 0x0 input");

  }
  if( D>2 && (!in3 || !in3_int) ){
    throw std::runtime_error( "unable to process 0x0 input");

  }

  // Check if output ndarray resides on the desired device
  //
  if( out ){
    if( alloc_device == CUNDA_CURRENT_DEVICE && out->get_device() != old_device ){
      out->set_device( old_device ); } // device copy
    else if( alloc_device == CUNDA_NDARRAY_DEVICE && out->get_device() != in1->get_device() ){
      out->set_device( in1->get_device() ); } // device copy
  }

  // Check if in_out ndarray resides on the desired device
  //
  if( out_idx > 0 && out_idx < 4 ){ 

   if( out_idx > D ){
      throw std::runtime_error( ">>>Internal error<<< :restore: array index out of range");

    }
   
   if( D == 1 ){
     if( in1->get_device() != in1_int->get_device() ){ 
       *in1 = *in1_int; } // device transfer by assignment
   }
   if( D == 2 ){
     if( out_idx == 1 && in1->get_device() != in1_int->get_device() ){ 
       *in1 = *in1_int; } // device transfer by assignment
     else if( out_idx == 2 && in2->get_device() != in2_int->get_device() ){ 
       *in2 = *in2_int; } // device transfer by assignment
   }
   if( D == 3 ){
     if( out_idx == 1 && in1->get_device() != in1_int->get_device() ){ 
       *in1 = *in1_int; } // device transfer by assignment
     else if( out_idx == 2 && in2->get_device() != in2_int->get_device() ){ 
       *in2 = *in2_int; } // device transfer by assignment
     else if( out_idx == 3 && in3->get_device() != in3_int->get_device() ){ 
       *in3 = *in3_int; } // device transfer by assignment
   }
  }
  else if( out_idx != 0 ){
    throw std::runtime_error( ">>>Internal error<<< :restore: illegal device specified");

  }

  // Check if internal array needs deletion (they do only if they were created in ::prepare()
  //
  if( in1->get_device() != in1_int->get_device() ){
    delete in1_int;
  }
  if( D>1 && in2->get_device() != in2_int->get_device() ){
    delete in2_int;
  }
  if( D>2 && in3->get_device() != in3_int->get_device() ){
    delete in3_int;
  }

  // Get current Cuda device
  int device;
  if( cudaGetDevice(&device) != cudaSuccess ) {
    throw cuda_error( "unable to get device no");

  }

  // Restore old device
  if( device != old_device && cudaSetDevice(old_device) != cudaSuccess) {
    throw cuda_error( "unable to restore device no");

  }
    

}

// Initialize static variables
//
static void initialize_static_variables()
{
  // This function is executed only once

  if( cudaGetDeviceCount( &num_devices ) != cudaSuccess) {
  	num_devices = 0;
    throw cuda_error( "Error: no Cuda devices present.");
  }

  int old_device;
  if( cudaGetDevice(&old_device) != cudaSuccess ) {
    throw std::runtime_error( "Error: unable to get device no");

  }

  warp_size = new int[num_devices];
  max_blockdim = new int[num_devices];
  max_griddim = new int[num_devices];
  handle = new cublasHandle_t[num_devices];

  if( !warp_size || !max_blockdim || !max_griddim || !handle ) {
    std::runtime_error( "Error: trivial malloc failed!"); // Do we really need to check this?? -DCHansen
  }

  for( int device=0; device<num_devices; device++ ){

    if( cudaSetDevice(device) != cudaSuccess ) {
      throw cuda_error( "Error: unable to set device no");

    }
    
    cudaDeviceProp deviceProp; 
    
    if( cudaGetDeviceProperties( &deviceProp, device ) != cudaSuccess) {
    	throw cuda_error("Error: unable to determine device properties.");

    }

    warp_size[device] = deviceProp.warpSize;
    max_blockdim[device] = deviceProp.maxThreadsDim[0];
    max_griddim[device] = deviceProp.maxGridSize[0];

    if (cublasCreate(&handle[device]) != CUBLAS_STATUS_SUCCESS) {
    	std::stringstream ss;
    	ss << "Error: unable to create cublas handle for device " << device << std::endl;
      throw cuda_error(ss.str());

    }

    cublasSetPointerMode( handle[device], CUBLAS_POINTER_MODE_HOST );

  }
  
  if( cudaSetDevice(old_device) != cudaSuccess ) {
    throw cuda_error( "Error: unable to restore device no");

  }

}

// Common block/grid configuration utility
//
static void setup_grid( unsigned int cur_device, unsigned int number_of_elements,
			dim3 *blockDim, dim3* gridDim, unsigned int num_batches=1 )
{
	initialize_static_variables();
  if( num_devices==0 ){
    throw cuda_error( "system device error");

  }

  // For small arrays we keep the block dimension fairly small
  *blockDim = dim3(256);
  *gridDim = dim3((number_of_elements+blockDim->x-1)/blockDim->x, num_batches);

  // Extend block/grid dimensions for large arrays
  if( gridDim->x > max_griddim[cur_device] ){
    blockDim->x = max_blockdim[cur_device];
    gridDim->x = (number_of_elements+blockDim->x-1)/blockDim->x;
  }

  if( gridDim->x > max_griddim[cur_device] ){
    gridDim->x = ((unsigned int)sqrt((float)number_of_elements)+blockDim->x-1)/blockDim->x;
    gridDim->y *= ((number_of_elements+blockDim->x*gridDim->x-1)/(blockDim->x*gridDim->x));
  }
   
  if( gridDim->x > max_griddim[cur_device] || gridDim->y > max_griddim[cur_device] ){

    throw cuda_error("Grid dimension larger than supported by device");
  }


}

// Common stride setup utility
//
template<class T> static void find_stride( cuNDArray<T> *in, unsigned int dim, 
					   unsigned int *stride, vector<unsigned int> *dims )
{
  *stride = 1;
  for( unsigned int i=0; i<in->get_number_of_dimensions(); i++ ){
    if( i != dim )
      dims->push_back(in->get_size(i));
    if( i < dim )
      *stride *= in->get_size(i);
  }
}

//
// Implementation of public utilities
//

// cAbs
//
template<class REAL, class T> __global__ void
cAbs_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    T val = in[idx];
    out[idx] = abs(val);
  }
}

// Abs
//
template<class T>
boost::shared_ptr< cuNDArray<typename realType<T>::type > >
abs( cuNDArray<T> *in,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
	typedef typename realType<T>::type REAL;
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Prepare 
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) cAbs_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// cNorm
//
template<class REAL, class T> __global__ void
cNorm_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    T val = in[idx];
    out[idx] = norm(val);
  }
}
/*
// cNorm
//
template<class REAL, class T>
boost::shared_ptr< cuNDArray<REAL> >
cNorm( cuNDArray<T> *in,
	    device alloc_device, device compute_device )
{
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Prepare
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    throw std::runtime_error( "cNorm: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get()); 
  if( out.get() != 0x0 ) cNorm_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
 
  // Restore 
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    throw std::runtime_error( "cNorm: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}

// Norm
//
template<class REAL, unsigned int D> __global__ void 
norm_kernel( typename reald<REAL,D>::Type *in, REAL *out,
		   unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    typename reald<REAL,D>::Type val = in[idx]; 
    out[idx] = norm<REAL,D>(val);
  }
}

// Norm
//
template<class REAL, unsigned int D>  
boost::shared_ptr< cuNDArray<REAL> > 
norm( cuNDArray<typename reald<REAL,D>::Type> *in,
	    device alloc_device, device compute_device )
{
  int cur_device, old_device; 
  cuNDArray< typename reald<REAL,D>::Type > *in_int;

  // Prepare
  prepare<1,typename reald<REAL,D>::Type,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    throw std::runtime_error( "norm: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    throw std::runtime_error( "norm: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) norm_kernel<REAL,D><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in_int->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,typename reald<REAL,D>::Type,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    throw std::runtime_error( "norm: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
 
  return out;
}

// Norm squared
//
template<class REAL, unsigned int D> __global__ 
void norm_squared_kernel( typename reald<REAL,D>::Type *in, REAL *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    typename reald<REAL,D>::Type val = in[idx]; 
    out[idx] = norm_squared<REAL,D>(val);
  }
} 

// Norm Squared
//
template<class REAL, unsigned int D>  
boost::shared_ptr< cuNDArray<REAL> > 
norm_squared( cuNDArray<typename reald<REAL,D>::Type> *in,
		    device alloc_device, device compute_device )
{
  int cur_device, old_device;
  cuNDArray< typename reald<REAL,D>::Type > *in_int;

  // Prepare
  prepare<1,typename reald<REAL,D>::Type,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    throw std::runtime_error( "norm_squared: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    throw std::runtime_error( "norm_squared: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) norm_squared_kernel<REAL,D><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  restore<1,typename reald<REAL,D>::Type,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    throw std::runtime_error( "norm_squared: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}
*/


// Sum
//
template<class T> __global__ void
sum_kernel( T *in, T *out,
		  unsigned int stride, 
		  unsigned int number_of_batches, 
		  unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
 
    T val = in[in_idx];
 
    for( unsigned int i=1; i<number_of_batches; i++ ) 
      val += in[i*stride+in_idx];

    out[idx] = val; 
  }
}

// Sum
//
template<class T>  
boost::shared_ptr< cuNDArray<T> > 
sum( cuNDArray<T> *in, unsigned int dim,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
  
  // Some validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    throw std::runtime_error("sum: underdimensioned.");
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    throw std::runtime_error( "sum: dimension out of range.");
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );
 
  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) sum_kernel<T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Expand
//
template<class T> __global__ void
expand_kernel( T *in, T *out,
		     unsigned int number_of_elements,
		     unsigned int new_dim_size )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = in[idx%(number_of_elements/new_dim_size)];
  }
}

// Expand
//
template<class T>  
boost::shared_ptr< cuNDArray<T> > 
expand( cuNDArray<T> *in, unsigned int new_dim_size,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
   
  unsigned int number_of_elements = in->get_number_of_elements()*new_dim_size;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );
 
  // Find element stride
  vector<unsigned int> dims = *in->get_dimensions(); 
  dims.push_back(new_dim_size);

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) expand_kernel<T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), number_of_elements, new_dim_size );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// SS
template<class REAL, class T> __inline__  __device__ REAL
_ss( unsigned int idx, T *in, unsigned int stride, unsigned int number_of_batches )
{
  unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
  REAL ss = REAL(0);

  for( unsigned int i=0; i<number_of_batches; i++ )
    ss += norm(in[i*stride+in_idx]);

  return ss;
}

// SS
template<class REAL, class T> __global__ void
ss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _ss<REAL,T>(idx, in, stride, number_of_batches);
  }
}

// SS
template<class REAL, class T>
boost::shared_ptr< cuNDArray<REAL> >
_ss( cuNDArray<T> *in, unsigned int dim,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    throw std::runtime_error( "ss: underdimensioned.");

  }

  if( dim > in->get_number_of_dimensions()-1 ){
    throw std::runtime_error( "ss: dimension out of range.");

  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if ( out.get() != 0x0 ) ss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
squaredNorm<float>( cuNDArray<float> *in, unsigned int dim,
		       cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _ss<float, float>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
squaredNorm<double>( cuNDArray<double> *in, unsigned int dim,
			 cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _ss<double, double>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
squaredNorm<float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
				      cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _ss<float, float_complext>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
squaredNorm<double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
		cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _ss<double, double_complext>(in, dim, alloc_device, compute_device);
}

// cSS
template<class REAL, class T> __global__ void
css_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    REAL ss = _ss<REAL,T>(idx, in, stride, number_of_batches); 

    out[idx].vec[0] = ss;
    out[idx].vec[1] = REAL(0);
  }
}

// cSS
template<class REAL>  
boost::shared_ptr< cuNDArray<complext<REAL> > >
css( cuNDArray<complext<REAL> > *in, unsigned int dim,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    std::runtime_error( "css: underdimensioned.");
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "css: dimension out of range." << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<complext<REAL> >( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(&dims);
  if ( out.get() != 0x0 ) css_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  restore<1,complext<REAL>,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );


  return out;
}

// RSS
template<class REAL, class T> __inline__  __device__ REAL
_rss( unsigned int idx, T *in, unsigned int stride, unsigned int number_of_batches )
{
  unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
  REAL rss = REAL(0);
  
  for( unsigned int i=0; i<number_of_batches; i++ ) 
    rss += norm(in[i*stride+in_idx]);
  
  rss = sqrt(rss); 

  return rss;
}

// RSS
template<class REAL, class T> __global__ void
rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _rss<REAL,T>(idx, in, stride, number_of_batches); 
  }
}

// RSS
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<REAL> > 
_rss( cuNDArray<T> *in, unsigned int dim,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    throw std::runtime_error( "rss: underdimensioned.");

  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    throw std::runtime_error( "rss: dimension out of range.");
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims); 
  if ( out.get() != 0x0 ) rss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}
/*
template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
norm<float>( cuNDArray<float> *in, unsigned int dim,
		cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _rss<float, float>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
norm<double>( cuNDArray<double> *in, unsigned int dim,
			  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _rss<double, double>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
norm<float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
				       cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _rss<float, float_complext>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
norm<double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
					 cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _rss<double, double_complext>(in, dim, alloc_device, compute_device);
}
*/
// cRSS
template<class REAL, class T> __global__ void
crss_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    REAL rss = _rss<REAL,T>(idx, in, stride, number_of_batches); 

    out[idx].vec[0] = rss;
    out[idx].vec[1] = REAL(0);
  }
}

// cRSS
template<class REAL>  
boost::shared_ptr< cuNDArray<complext<REAL> > >
crss( cuNDArray<complext<REAL> > *in, unsigned int dim,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "crss: underdimensioned." << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "crss: dimension out of range." << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ;

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<complext<REAL> >( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(&dims);
  if ( out.get() != 0x0 ) crss_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  restore<1,complext<REAL>,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );
  
  return out;
}


/*
// reciprocal RSS
template<class REAL, class T> __global__ void
reciprocal_rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = 1/(_rss<REAL,T>(idx, in, stride, number_of_batches));
  }
}

// Reciprocal RSS
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<REAL> > 
_reciprocal_rss( cuNDArray<T> *in, unsigned int dim,
		       device alloc_device, device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "reciprocal_rss: underdimensioned." << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "reciprocal_rss: dimension out of range." << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if ( out.get() != 0x0 ) reciprocal_rss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );
  
  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
reciprocal_rss<float,float>( cuNDArray<float> *in, unsigned int dim,
				   device alloc_device, device compute_device )
{
  return _reciprocal_rss<float, float>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
reciprocal_rss<double,double>( cuNDArray<double> *in, unsigned int dim,
				     device alloc_device, device compute_device )
{
  return _reciprocal_rss<double, double>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
reciprocal_rss<float,float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
						  device alloc_device, device compute_device )
{
  return _reciprocal_rss<float, float_complext>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
reciprocal_rss<double,double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
						    device alloc_device, device compute_device )
{
  return _reciprocal_rss<double, double_complext>(in, dim, alloc_device, compute_device);
}

// cReciprocal RSS
template<class REAL, class T> __global__ void
creciprocal_rss_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    REAL reciprocal_rss = 1/(_rss<REAL,T>(idx, in, stride, number_of_batches));

    out[idx].vec[0] = reciprocal_rss;
    out[idx].vec[1] = REAL(0);
  }
}

// cReciprocal RSS
template<class REAL>  
boost::shared_ptr< cuNDArray<complext<REAL> > >
creciprocal_rss( cuNDArray<complext<REAL> > *in, unsigned int dim,
		       device alloc_device, device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    throw std::runtime_error( "creciprocal_rss: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "creciprocal_rss: underdimensioned." << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "creciprocal_rss: dimension out of range." << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    throw std::runtime_error( "creciprocal_rss: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<complext<REAL> >( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(&dims);
  if ( out.get() != 0x0 ) creciprocal_rss_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  restore<1,complext<REAL>,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    throw std::runtime_error( "creciprocal_rss: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float_complext> >
reciprocal_rss<float_complext, float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
								  device alloc_device, device compute_device )
{
  return creciprocal_rss<float>( in, dim, alloc_device, compute_device );
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double_complext> >
reciprocal_rss<double_complext, double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
								    device alloc_device, device compute_device )
{
  return creciprocal_rss<double>( in, dim, alloc_device, compute_device );
}
*/

// Build correlation matrix
template<class REAL, class T> __global__ void
correlation_kernel( T *in, T *corrm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.y;

  if( p < num_elements ){
    for( unsigned int j=0; j<i; j++){
      T tmp = in[i*num_elements+p]*conj(in[j*num_elements+p]);
      corrm[(j*num_batches+i)*num_elements+p] = tmp;
      corrm[(i*num_batches+j)*num_elements+p] = conj(tmp);
    }
    T tmp = in[i*num_elements+p];
    corrm[(i*num_batches+i)*num_elements+p] = tmp*conj(tmp);
  }
}

// Build correlation matrix
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<T> >
_correlation( cuNDArray<T> *in,
		    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "correlation: underdimensioned." << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }
 
  unsigned int number_of_batches = in->get_size(in->get_number_of_dimensions()-1);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  dim3 blockDim(((max_blockdim[old_device]/number_of_batches)/warp_size[old_device])*warp_size[old_device], number_of_batches);

  if( blockDim.x == 0 ){
    cout << endl << "correlation: correlation dimension exceeds device capacity." << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }
  
  dim3 gridDim((number_of_elements+blockDim.x-1)/blockDim.x);

  // Invoke kernel
  vector<unsigned int> dims = *in->get_dimensions(); dims.push_back(number_of_batches);
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) correlation_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
correlation<float>( cuNDArray<float> *data,
			  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _correlation<float,float>(data, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> >
correlation<double>( cuNDArray<double> *data,
			   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _correlation<double,double>(data, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float_complext> >
correlation<float_complext>( cuNDArray<float_complext> *data,
					 cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _correlation<float,float_complext>(data, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double_complext> >
correlation<double_complext>( cuNDArray<double_complext> *data,
					  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _correlation<double,double_complext>(data, alloc_device, compute_device);
}

// Real to complext
template<class REAL> __global__ void
real_to_complext_kernel( REAL *in, complext<REAL> *out, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < num_elements ){
    complext<REAL> z;
    z.vec[0] = in[idx];
    z.vec[1] = REAL(0);
    out[idx] = z;
  }
}

// Convert real to complext
template<class REAL>  
boost::shared_ptr< cuNDArray<complext<REAL> > >
real_to_complext( cuNDArray<REAL> *in,
			cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
 
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) real_to_complext_kernel<REAL><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements());
  
  CHECK_FOR_CUDA_ERROR();
  
  // Restore
  restore<1,REAL,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// complext to real by cropping the imaginary component 
template<class REAL> __global__ void
complext_to_real_kernel( complext<REAL> *in, REAL *out, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < num_elements ){
    out[idx] = in[idx].vec[0];
  }
}

// Convert complext to real by cropping the imaginary component 
template<class REAL>  
boost::shared_ptr< cuNDArray<REAL> > 
complext_to_real( cuNDArray<complext<REAL> > *in,
			cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ;
 
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());  
  if( out.get() != 0x0 ) complext_to_real_kernel<REAL><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements());
  
  CHECK_FOR_CUDA_ERROR();
  
  // Restore
  restore<1,complext<REAL>,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Downsample
template<class REAL, unsigned int D> __global__ void
downsample_kernel( REAL *in, REAL *out,
			 typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
			 unsigned int num_elements, unsigned int num_batches )
{
  // We have started a thread for each output element
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  const unsigned int frame_offset = idx/num_elements;
  
  if( idx < num_elements*num_batches ){

    const typename uintd<D>::Type co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
    const typename uintd<D>::Type co_in = co_out << 1;

    const typename uintd<D>::Type twos = to_vector_td<unsigned int,D>(2);
    const unsigned int num_adds = 1 << D;
    unsigned int actual_adds = 0;

    REAL res = REAL(0);

    for( unsigned int i=0; i<num_adds; i++ ){
      const typename uintd<D>::Type local_co = idx_to_co<D>( i, twos );
      if( weak_greater_equal( local_co, matrix_size_out ) ) continue; // To allow array dimensions of 1
      const unsigned int in_idx = co_to_idx<D>(co_in+local_co, matrix_size_in)+frame_offset*prod(matrix_size_in); 
      actual_adds++;
      res += in[in_idx];
    }
    
    out[idx] = res/REAL(actual_adds);
  }
}

// Downsample
template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<REAL> > 
downsample( cuNDArray<REAL> *in,
		  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    throw std::runtime_error( "downsample: the number of array dimensions should be at least D");

  }
  
  for( unsigned int d=0; d<D; d++ ){
    if( (in->get_size(d)%2) == 1 && in->get_size(d) != 1 ){
      throw std::runtime_error( "downsample: uneven array dimensions larger than one not accepted");
    }
  }
  
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = matrix_size_in >> 1;

  for( unsigned int d=0; d<D; d++ ){
    if( matrix_size_out[d] == 0 ) 
      matrix_size_out[d] = 1;
  }
  
  unsigned int number_of_elements = prod(matrix_size_out);
  unsigned int number_of_batches = 1;

  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;

  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches ) ;
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out); 
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    downsample_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Nearest neighbor upsampling
template<class REAL, unsigned int D> __global__ void
upsample_nn_kernel( REAL *in, REAL *out,
		       typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
		       unsigned int num_elements, unsigned int num_batches )
{
  // We have started a thread for each output element
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx < num_elements*num_batches ){    
    const unsigned int frame_idx = idx/num_elements;
    const typename uintd<D>::Type co_out = idx_to_co<D>( idx-frame_idx*num_elements, matrix_size_out );
    const typename uintd<D>::Type co_in = co_out >> 1;
    out[idx] = in[co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in)];
  }
}

// Nearest neighbor upsampling
template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<REAL> > 
upsample_nn( cuNDArray<REAL> *in,
		   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    throw std::runtime_error( "upsample: the number of array dimensions should be at least D" );

  }
    
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = matrix_size_in << 1;

  unsigned int number_of_elements = prod(matrix_size_out);
  unsigned int number_of_batches = 1;

  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;

  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches );
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out); 
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    upsample_nn_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Utility to check if all neighbors required for the linear interpolation exists
// ... do not include dimensions of size 1

template<class REAL, unsigned int D> __device__ 
bool is_border_pixel( typename uintd<D>::Type co, typename uintd<D>::Type dims )
{
  for( unsigned int dim=0; dim<D; dim++ ){
    if( dims[dim] > 1 && ( co[dim] == 0 || co[dim] == (dims[dim]-1) ) )
      return true;
  }
  return false;
}

// Linear upsampling
template<class REAL, unsigned int D> __global__ void
upsample_lin_kernel( REAL *in, REAL *out,
		       typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
		       unsigned int num_elements, unsigned int num_batches )
{
  // We have started a thread for each output element
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx < num_elements*num_batches ){

    REAL res = REAL(0);

    const unsigned int num_neighbors = 1 << D;
    const unsigned int frame_idx = idx/num_elements;
    const typename uintd<D>::Type co_out = idx_to_co<D>( idx-frame_idx*num_elements, matrix_size_out );

    // We will only proceed if all neighbours exist (this adds a zero-boundary to the upsampled image/vector field)
    //
    
    if( !is_border_pixel<REAL,D>(co_out, matrix_size_out) ){
      
      for( unsigned int i=0; i<num_neighbors; i++ ){
	
	// Determine coordinate of neighbor in input
	//

	const typename uintd<D>::Type twos = to_vector_td<unsigned int,D>(2);
	const typename uintd<D>::Type stride = idx_to_co<D>( i, twos );

	if( weak_greater_equal( stride, matrix_size_out ) ) continue; // To allow array dimensions of 1

	// Be careful about dimensions of size 1
	typename uintd<D>::Type ones = to_vector_td<unsigned int,D>(1);
	for( unsigned int d=0; d<D; d++ ){
	  if( matrix_size_out[d] == 1 )
	    ones[d] = 0;
	}
	typename uintd<D>::Type co_in = ((co_out-ones)>>1)+stride;
	
	// Read corresponding pixel value
	//
	
	const unsigned int in_idx = co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in); 
	REAL value = in[in_idx];
	
	// Determine weight
	//
	
	REAL weight = REAL(1);
	
	for( unsigned int dim=0; dim<D; dim++ ){	  
	  if( matrix_size_in[dim] > 1 ){
	    if( stride.vec[dim] == (co_out.vec[dim]%2) ) {
	      weight *= REAL(0.25);
	    }
	    else{
	      weight *= REAL(0.75);
	    }
	  }
	}
	
	// Accumulate result
	//
	
	res += weight*value;
      }
    }
    out[idx] = res;
  }
}

// Linear interpolation upsampling
template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<REAL> > 
upsample_lin( cuNDArray<REAL> *in,
		    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    throw std::runtime_error( "upsample: the number of array dimensions should be at least D");
  }
    
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = matrix_size_in << 1;

  for( unsigned int d=0; d<D; d++ ){
    if( matrix_size_in[d] == 1 )
      matrix_size_out[d] = 1;
  }
  
  unsigned int number_of_elements = prod(matrix_size_out);
  unsigned int number_of_batches = 1;

  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;

  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches );
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out); 
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    upsample_lin_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}




// threshold
template<class T> __global__
void threshold_min_kernel( T min, T value, T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    if( real(in_out[idx]) < real(min) ) 
      in_out[idx] = value;
  }
}

//Threshold
template<class T>
void threshold_min( T min, cuNDArray<T> *in_out, T value, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  threshold_min_kernel<T><<< gridDim, blockDim >>>(min, value, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}

// threshold
template<class T> __global__
void threshold_max_kernel( T max, T value, T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    if( real(in_out[idx]) > real(max) ) 
      in_out[idx] = value;
  }
}


//Threshold
template<class T>
void threshold_max( T max, cuNDArray<T> *in_out, T value, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  threshold_max_kernel<T><<< gridDim, blockDim >>>(max, value, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}

// threshold
template<class T> __global__
void threshold_min_kernel2( T * min, T value, T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
   
  if( idx<number_of_elements ){
    if( real(in_out[idx]) < real(min[idx]) ) 
      in_out[idx] = value;
  }
}

//Threshold
template<class T>
void threshold_min( cuNDArray<T> * min, cuNDArray<T> *in_out, T value, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;
  cuNDArray<T> *min_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int, min, &min_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  threshold_min_kernel2<T><<< gridDim, blockDim >>>(min_int->get_data_ptr(), value, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in_out, in_out_int, 1,compute_device,0x0,min,min_int );


}

// threshold
template<class T> __global__
void threshold_max_kernel2( T * min, T value, T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
    if( real(in_out[idx]) > real(min[idx]) ) 
      in_out[idx] = value;
  }
}

//Threshold
template<class T>
void threshold_max( cuNDArray<T> * max, cuNDArray<T> *in_out, T value, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;
  cuNDArray<T> *max_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int, max, &max_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  threshold_max_kernel2<T><<< gridDim, blockDim >>>(max_int->get_data_ptr(), value, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in_out, in_out_int, 1,compute_device,0,max,max_int );


}


// threshold
template<class T> __global__
void threshold_amin_kernel( T * min, T value, T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
   
  if( idx<number_of_elements ){
    if( abs(in_out[idx]) < real(min[idx]) ) 
      in_out[idx] = value;
  }
}

//Threshold
template<class T>
void threshold_amin( cuNDArray<T> * min, cuNDArray<T> *in_out, T value, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;
  cuNDArray<T> *min_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int, min, &min_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  threshold_amin_kernel<T><<< gridDim, blockDim >>>(min_int->get_data_ptr(), value, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in_out, in_out_int, 1,compute_device,0,min,min_int );


}

// Reciprocal square root
template<class T> __global__ 
void reciprocal_sqrt_kernel( T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = rsqrt(in_out[idx]);
  }
}

// Square root
template<class T> 
void reciprocal_sqrt( cuNDArray<T> *in_out,
			    cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  reciprocal_sqrt_kernel<T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), in_out->get_number_of_elements() );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );

}

// Normalized RSS
template<class REAL, class T> __global__ void
rss_normalize_kernel( T *in_out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    REAL reciprocal_rss = 1/(_rss<REAL,T>(idx, in_out, stride, number_of_batches));
 
    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

    for( unsigned int i=0; i<number_of_batches; i++ ) {
      T out = in_out[i*stride+in_idx];
      out *= reciprocal_rss; // complex-scalar multiplication (element-wise operator)
      in_out[i*stride+in_idx] = out; 
    } 
  }
}

// Normalized RSS
template<class T>
void rss_normalize( cuNDArray<T> *in_out, unsigned int dim,
			  cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int );

  // Validity checks
  if( !(in_out->get_number_of_dimensions()>1) ){
    throw std::runtime_error( "rss_normalize: underdimensioned.");

  }
 
  if( dim > in_out->get_number_of_dimensions()-1 ){
  	throw std::runtime_error("rss_normalize: dimension out of range.");

  }

  unsigned int number_of_batches = in_out->get_size(dim);
  unsigned int number_of_elements = in_out->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );
  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in_out, dim, &stride, &dims );

  // Invoke kernel
  rss_normalize_kernel<typename realType<T>::type,T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}



// Scale conjugate w. non conjugate
template<class S, class T> __global__ 
void scale_conj_kernel( S *a, T *x, unsigned int number_of_elements, unsigned int number_of_batches )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements*number_of_batches ){
    unsigned int frame_offset = idx/number_of_elements;
    S in_a = a[idx-frame_offset*number_of_elements];
    T in_x = x[idx];
    x[idx] = conj(in_a)*in_x;
  }
}

// Scale conjugate w. non conjugate
template<class T> 
bool scale_conj( cuNDArray<T> *a, cuNDArray<T> *x,
		  cuNDA_device compute_device )
{
  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch, cannot scale" << endl;
    return false;
  }
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *a_int, *x_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, a, &a_int, x, &x_int );

  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, num_batches );
  
  // Invoke kernel
  scale_conj_kernel<T,T><<< gridDim, blockDim >>> ( a_int->get_data_ptr(), x_int->get_data_ptr(), number_of_elements, num_batches );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, a, a_int, 2, compute_device, 0x0, x, x_int );

  return true;
}



// Normalize (float)
template<> EXPORTGPUCORE
float normalize<float>( cuNDArray<float> *data, float new_max, cuNDA_device compute_device )
{
  initialize_static_variables();


  unsigned int number_of_elements = data->get_number_of_elements();

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<float> *data_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,float,dummy,dummy>( compute_device, &cur_device, &old_device, data, &data_int );

  // Find the maximum value in the array
  int max_idx;
  cublasIsamax( handle[cur_device], number_of_elements, data_int->get_data_ptr(), 1, &max_idx );
  cudaThreadSynchronize();
  
  // Copy that value back to host memory
  float max_val;
  cudaMemcpy(&max_val, (data_int->get_data_ptr()+max_idx-1), sizeof(float), cudaMemcpyDeviceToHost);

  // Scale the array
  float scale = abs(new_max/max_val);
  cublasSscal( handle[cur_device], number_of_elements, &scale, data_int->get_data_ptr(), 1 );

  // Restore
  restore<1,float,dummy,dummy,dummy>( old_device, data, data_int, 1, compute_device );

  CHECK_FOR_CUDA_ERROR();
  return scale;
}

// Normalize (double)
template<> EXPORTGPUCORE
double normalize<double>( cuNDArray<double> *data, double new_max, cuNDA_device compute_device )
{
  initialize_static_variables();

  unsigned int number_of_elements = data->get_number_of_elements();

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<double> *data_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,double,dummy,dummy>( compute_device, &cur_device, &old_device, data, &data_int );

  // Find the maximum value in the array
  int max_idx;
  cublasIdamax( handle[cur_device], number_of_elements, data_int->get_data_ptr(), 1, &max_idx );
  cudaThreadSynchronize();

  // Copy that value back to host memory
  double max_val;
  cudaMemcpy(&max_val, (data_int->get_data_ptr()+max_idx-1), sizeof(double), cudaMemcpyDeviceToHost);

  // Scale the array
  double scale = abs(new_max/max_val);
  cublasDscal( handle[cur_device], number_of_elements, &scale, data_int->get_data_ptr(), 1 );
  
  // Restore
  restore<1,double,dummy,dummy,dummy>( old_device, data, data_int, 1, compute_device );
  
  CHECK_FOR_CUDA_ERROR();
  return scale;
}


// Crop
template<class T, unsigned int D> __global__ void
crop_kernel( typename uintd<D>::Type offset, typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out,
		   T *in, T *out, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  const unsigned int frame_offset = idx/num_elements;

  if( idx < num_elements*num_batches ){
    const typename uintd<D>::Type co = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
    const typename uintd<D>::Type co_os = offset + co;
    const unsigned int in_idx = co_to_idx<D>(co_os, matrix_size_in)+frame_offset*prod(matrix_size_in); 
    out[idx] = in[in_idx];
  }
}

// Crop
template<class T, unsigned int D> EXPORTGPUCORE
void crop( typename uintd<D>::Type offset,
	    cuNDArray<T> *in, cuNDArray<T> *out,
	    cuNDA_device compute_device )
{
  if( in == 0x0 || out == 0x0 ){

    throw std::runtime_error("crop: 0x0 ndarray provided");
  }

  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
  	throw std::runtime_error("crop: image dimensions mismatch");

  }

  if( in->get_number_of_dimensions() < D ){
    std::stringstream ss;
    ss << "crop: number of image dimensions should be at least " << D;

    throw std::runtime_error(ss.str());
  }

  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *out->get_dimensions() );
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }

  if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
    throw std::runtime_error( "crop: cropping size mismatch");

  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int );
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

  // Invoke kernel
  crop_kernel<T,D><<< gridDim, blockDim >>>
    ( offset, matrix_size_in, matrix_size_out, in_int->get_data_ptr(), out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int );

}

// Expand and zero fill
template<class T, unsigned int D> __global__ void
expand_with_zero_fill_kernel( typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out,
				    T *in, T *out, unsigned int number_of_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  const unsigned int frame_offset = idx/num_elements;

  if( idx < num_elements*number_of_batches ){

    const typename uintd<D>::Type co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
    const typename uintd<D>::Type offset = (matrix_size_out-matrix_size_in)>>1;
    T _out;
    bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));

    if( inside )
      _out = in[co_to_idx<D>(co_out-offset, matrix_size_in)+frame_offset*prod(matrix_size_in)];
    else{      
      _out = T(0);
    }

    out[idx] = _out;
  }
}

// Expand and zero fill
template<class T, unsigned int D> 
void expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out,
				  cuNDA_device compute_device )
{ 
  if( in == 0x0 || out == 0x0 ){
  	throw std::runtime_error("zero_fill: 0x0 ndarray provided");

  }

  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
  	throw std::runtime_error("zero_fill: image dimensions mismatch");

  }

  if( in->get_number_of_dimensions() < D ){
  	std::stringstream ss;
    ss << "zero_fill: number of image dimensions should be at least " << D;
    throw std::runtime_error(ss.str());

  }

  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *out->get_dimensions() );
  
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }

  if( weak_greater(matrix_size_in,matrix_size_out) ){
    std::runtime_error("expand: size mismatch, cannot expand");

  }
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );
 
  // Invoke kernel
  expand_with_zero_fill_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, in_int->get_data_ptr(), out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int );

}

// Zero fill border (rectangular)
template<class T, unsigned int D> __global__ void
zero_fill_border_kernel( typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out,
			       T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
    const typename uintd<D>::Type offset = (matrix_size_out-matrix_size_in)>>1;
    if( weak_less( co_out, offset ) || weak_greater_equal( co_out, matrix_size_in+offset ) ){
      T zero = T(0);
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
void zero_fill_border( typename uintd<D>::Type matrix_size_in, cuNDArray<T> *in_out,
			     cuNDA_device compute_device )
{ 
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *in_out->get_dimensions() );
 
  if( weak_greater(matrix_size_in, matrix_size_out) ){
    std::runtime_error( "zero_fill: size mismatch, cannot zero fill");

  }
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
    number_of_batches *= in_out->get_size(d);
  }

 // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim );
 
  // Invoke kernel
  zero_fill_border_kernel<T,D><<< gridDim, blockDim >>>
  		( matrix_size_in, matrix_size_out, in_out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}

// Zero fill border (circular)
template<class REAL, class T, unsigned int D> __global__ void
zero_fill_border_kernel( REAL radius, typename intd<D>::Type dims, T *image,
			       unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  const unsigned int frame_offset = idx/num_elements;
  
  if( idx < num_elements*num_batches ){
    const typename intd<D>::Type co = idx_to_co<D>( idx-frame_offset*num_elements, dims ) - (dims>>1);
    if( REAL(norm_squared(co)) > radius*radius )
      image[idx] = T(0);
  }
}

// Zero fill border (circular, 2D)
template<class REAL, class T, unsigned int D> 
void zero_fill_border( REAL radius, cuNDArray<T> *in_out,
			     cuNDA_device compute_device )
{

  
  unsigned int number_of_batches = 1;
  for( unsigned int d=2; d<in_out->get_number_of_dimensions(); d++ ){
    number_of_batches *= in_out->get_size(d);
  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int );

  typename intd<D>::Type dims = vector_to_intd<D>(*in_out->get_dimensions());

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim, number_of_batches );
 
  // Invoke kernel
  zero_fill_border_kernel<REAL,T,D><<< gridDim, blockDim >>>
    ( radius, dims, in_out_int->get_data_ptr(), number_of_batches, prod(dims) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}

// Shrinkage
//

template<class REAL, class T> __global__ void 
shrink1_kernel( REAL gamma, T *in, T *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T in_val = in[idx]; 
    REAL in_norm = abs<REAL>(in_val);
    T _res;
    if( in_norm > REAL(0) )
      _res =  in_val/in_norm;
    else
      _res = T(0);
    REAL maximum = max( in_norm-gamma, REAL(0) );
    T res = maximum*_res;

    out[idx] = res;
  }
}

template<class REAL, class T> EXPORTGPUCORE
void shrink1( REAL gamma, cuNDArray<T> *in, cuNDArray<T> *out )
{
  // TODO: multi-device handling

  if( !in || !out ){
    throw std::runtime_error( "shrink1: 0x0 arrays not accepted" );

  }

  if( in->get_number_of_elements() != out->get_number_of_elements() ){
    throw std::runtime_error( "shrink1: i/o arrays must have an identical number of elements");
  }
  
  // Get current Cuda device
  int cur_device;
  if( cudaGetDevice(&cur_device) != cudaSuccess ) {
    throw std::runtime_error( "shrink1 : unable to get device no");
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );
  
  // Invoke kernel
  shrink1_kernel<REAL,T><<< gridDim, blockDim >>>( gamma, in->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
  

}

template<class REAL, class T> __global__ void 
shrinkd_kernel( REAL gamma, REAL *s_k, T *in, T *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T in_val = in[idx]; 
    REAL s_k_val = s_k[idx];
    T _res;
    if( s_k_val > 0 )
      _res =  in_val/s_k_val;
    else
      _res = T(0);
    REAL maximum = max( s_k_val-gamma, REAL(0) );
    T res = maximum*_res;

    out[idx] = res;
  }
}

template<class REAL, class T> EXPORTGPUCORE
void shrinkd( REAL gamma, cuNDArray<REAL> *s_k, cuNDArray<T> *in, cuNDArray<T> *out )
{
  // TODO: multi-device handling

  if( !in || !out || !s_k ){
    throw std::runtime_error( "shrinkd: 0x0 arrays not accepted");

  }

  if( in->get_number_of_elements() != out->get_number_of_elements() ){
    throw std::runtime_error( "shrinkd: i/o arrays must have an identical number of elements");

  }

  if( in->get_number_of_elements() != s_k->get_number_of_elements() ){
    throw std::runtime_error( "shrinkd: i/o arrays must have an identical number of elements");

  }
  
  // Get current Cuda device
  int cur_device;
  if( cudaGetDevice(&cur_device) != cudaSuccess ) {
    throw std::runtime_error( "shrinkd : unable to get device no");

  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );
  
  // Invoke kernel
  shrinkd_kernel<REAL,T><<< gridDim, blockDim >>>( gamma, s_k->get_data_ptr(), in->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();

}

// Mirror, but keep the origin unchanged
template<class T, unsigned int D> __global__ void
origin_mirror_kernel( typename uintd<D>::Type matrix_size, typename uintd<D>::Type origin, T *in, T *out, bool zero_fill )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < prod(matrix_size) ){

    typename uintd<D>::Type in_co = idx_to_co<D>( idx, matrix_size );
    typename uintd<D>::Type out_co = matrix_size-in_co;
    
    bool wrap = false;
    for( unsigned int d=0; d<D; d++ ){
      if( out_co.vec[d] == matrix_size.vec[d] ){
	out_co.vec[d] = 0;
	wrap = true;
      }
    }
    
    const unsigned int in_idx = co_to_idx<D>(in_co, matrix_size);
    const unsigned int out_idx = co_to_idx<D>(out_co, matrix_size);

    if( wrap && zero_fill )
      out[out_idx] = T(0);
    else
      out[out_idx] = in[in_idx];
  }
}

// Mirror around the origin -- !! leaving the origin unchanged !!
// This creates empty space "on the left" that can be filled by zero (default) or the left-over entry.
template<class T, unsigned int D> EXPORTGPUCORE
void origin_mirror( cuNDArray<T> *in, cuNDArray<T> *out, bool zero_fill, cuNDA_device compute_device )
{
  if( in == 0x0 || out == 0x0 ){
  	throw std::runtime_error( "origin_mirror: 0x0 ndarray provided");

  }

  if( !in->dimensions_equal(out) ){
  	throw std::runtime_error("origin_mirror: image dimensions mismatch");

  }
  
  if( in->get_number_of_dimensions() != D ){
  	std::stringstream ss;
    ss << "origin_mirror: number of image dimensions is not " << D;
    throw std::runtime_error(ss.str());
  }

  typename uintd<D>::Type matrix_size = vector_to_uintd<D>( *in->get_dimensions() );
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int );
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size), &blockDim, &gridDim );

  // Invoke kernel
  origin_mirror_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size, matrix_size>>1, in_int->get_data_ptr(), out_int->get_data_ptr(), zero_fill );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int ) ;

}



// Minimum
//
template<class T> __global__ 
void minimum_kernel( T* in1,T* in2, T* out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    out[idx]=min(in1[idx],in2[idx]);
  }
} 


// Minimum
//
template<class T>  
boost::shared_ptr< cuNDArray<T> > 
minimum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device;
  cuNDArray<T> *in1_int;
  cuNDArray<T> *in2_int;


  if ( in1->get_number_of_elements() !=  in2->get_number_of_elements()){
    throw std::runtime_error( "minimum: input arrays have different number of elements");

  }
  // Prepare 
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in1, &in1_int, in2, &in2_int);

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in1->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(in1->get_dimensions().get());
  if( out.get() != 0x0 ) minimum_kernel<T><<< gridDim, blockDim >>>( in1_int->get_data_ptr(), in2_int->get_data_ptr(),out->get_data_ptr(), in1->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore 
  // Restore
  restore<2,T,T,T,dummy>( old_device, in1, in1_int, 0, compute_device, out.get(), in2, in2_int );

  return out;
}



// Maximum
//
template<class T> __global__ 
void maximum_kernel( T* in1,T* in2, T* out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    out[idx]=max(in1[idx],in2[idx]);
  }
} 


// Minimum
//
template<class T>  
boost::shared_ptr< cuNDArray<T> > 
maximum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device;
  cuNDArray<T> *in1_int;
  cuNDArray<T> *in2_int;


  if ( in1->get_number_of_elements() !=  in2->get_number_of_elements()){
    throw std::runtime_error( "maximum: input arrays have different number of elements");
  }
  // Prepare 
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in1, &in1_int, in2, &in2_int  ) ;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in1->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(in1->get_dimensions().get());
  if( out.get() != 0x0 ) maximum_kernel<T><<< gridDim, blockDim >>>( in1_int->get_data_ptr(), in2_int->get_data_ptr(),out->get_data_ptr(), in1->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore 
  // Restore
  restore<2,T,T,T,dummy>( old_device, in1, in1_int, 0, compute_device, out.get(), in2, in2_int );

  return out;
}


//
// Instantiation
//

// A few functions have integer support

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<int> > 
sum<int>( cuNDArray<int>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<1>::Type> > 
sum<intd<1>::Type >( cuNDArray<intd<1>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<2>::Type> > 
sum<intd<2>::Type >( cuNDArray<intd<2>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<3>::Type> > 
sum<intd<3>::Type >( cuNDArray<intd<3>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<4>::Type> > 
sum<intd<4>::Type >( cuNDArray<intd<4>::Type >*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<unsigned int> > 
sum<unsigned int>( cuNDArray<unsigned int>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<1>::Type> > 
sum<uintd<1>::Type>( cuNDArray<uintd<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<2>::Type> > 
sum<uintd<2>::Type>( cuNDArray<uintd<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<3>::Type> > 
sum<uintd<3>::Type>( cuNDArray<uintd<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<4>::Type> > 
sum<uintd<4>::Type>( cuNDArray<uintd<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );


template EXPORTGPUCORE void
crop<int,1>( uintd1, cuNDArray<int>*, cuNDArray<int>*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<int,1>,1>( uintd1, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,2>,1>( uintd1, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,3>,1>( uintd1, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,4>,1>( uintd1, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<int,1>,2>( uintd2, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,2>,2>( uintd2, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,3>,2>( uintd2, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,4>,2>( uintd2, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<int,1>,3>( uintd3, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,2>,3>( uintd3, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,3>,3>( uintd3, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,4>,3>( uintd3, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<int,1>,4>( uintd4, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,2>,4>( uintd4, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,3>,4>( uintd4, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<int,4>,4>( uintd4, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<unsigned int,1>( uintd1, cuNDArray<unsigned int>*, cuNDArray<unsigned int>*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<unsigned int,1>,1>( uintd1, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,2>,1>( uintd1, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,3>,1>( uintd1, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,4>,1>( uintd1, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<unsigned int,1>,2>( uintd2, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,2>,2>( uintd2, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,3>,2>( uintd2, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,4>,2>( uintd2, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<unsigned int,1>,3>( uintd3, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,2>,3>( uintd3, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,3>,3>( uintd3, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,4>,3>( uintd3, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<unsigned int,1>,4>( uintd4, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,2>,4>( uintd4, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,3>,4>( uintd4, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<unsigned int,4>,4>( uintd4, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

// Instanciation -- single precision

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
sum<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<complext<float> > >
sum<complext<float> >( cuNDArray<complext<float> >*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<1>::Type> > 
sum<floatd<1>::Type>( cuNDArray<floatd<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<2>::Type> > 
sum<floatd<2>::Type>( cuNDArray<floatd<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<3>::Type> >
 sum<floatd<3>::Type>( cuNDArray<floatd<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<4>::Type> >
 sum<floatd<4>::Type>( cuNDArray<floatd<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
expand<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
expand<float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device);


template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
abs<float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
abs<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
minimum<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
maximum<float>( cuNDArray<float>*, cuNDArray<float>*,cuNDA_device, cuNDA_device );
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm<float,1>( cuNDArray<floatd<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm<float,2>( cuNDArray<floatd<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm<float,3>( cuNDArray<floatd<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm<float,4>( cuNDArray<floatd<4>::Type>*, cuNDA_device, cuNDA_device );
*/

/*template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
norm_squared<float,float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
 norm_squared<float,float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );*/
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm_squared<float,1>( cuNDArray<floatd<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm_squared<float,2>( cuNDArray<floatd<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm_squared<float,3>( cuNDArray<floatd<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
norm_squared<float,4>( cuNDArray<floatd<4>::Type>*, cuNDA_device, cuNDA_device );
*/
/*
template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
ss<float,float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
ss<float,float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
ss<float_complext, float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
/*
template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
rss<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
rss<float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
/*
template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
rss<float_complext, float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
/*
template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
reciprocal_rss<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
reciprocal_rss<float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
/*
template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
reciprocal_rss<float_complext, float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
correlation<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
correlation<float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );

//template EXPORTGPUCORE void axpy<float>( cuNDArray<float>*, cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void
crop<float,1>( uintd1, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
crop<float,2>( uintd2, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );

template EXPORTGPUCORE void
crop<complext<float>,1>( uintd1, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<complext<float>,2>( uintd2, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<complext<float>,3>( uintd3, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<complext<float>,4>( uintd4, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<float,1>,1>( uintd1, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,2>,1>( uintd1, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,3>,1>( uintd1, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,4>,1>( uintd1, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<float,1>,2>( uintd2, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,2>,2>( uintd2, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,3>,2>( uintd2, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,4>,2>( uintd2, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<float,1>,3>( uintd3, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,2>,3>( uintd3, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,3>,3>( uintd3, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,4>,3>( uintd3, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<float,1>,4>( uintd4, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,2>,4>( uintd4, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,3>,4>( uintd4, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<float,4>,4>( uintd4, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );


template EXPORTGPUCORE void
expand_with_zero_fill<float,1>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<float_complext,1>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void
expand_with_zero_fill<float,2>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<float_complext,2>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void
expand_with_zero_fill<float,3>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<float_complext,3>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void
expand_with_zero_fill<float,4>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<float_complext,4>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
real_to_complext<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
complext_to_real<float>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
downsample<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
downsample<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
downsample<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
downsample<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_nn<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_nn<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_nn<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_nn<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_lin<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_lin<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_lin<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
upsample_lin<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

//template EXPORTGPUCORE void clear<float>( cuNDArray<float>*,float, cuNDA_device );
//template EXPORTGPUCORE void clear<float_complext>( cuNDArray<float_complext>*,float_complext, cuNDA_device );

//template EXPORTGPUCORE void reciprocal<float>( cuNDArray<float>*, cuNDA_device );
//template EXPORTGPUCORE void reciprocal<float_complext>( cuNDArray<float_complext>*, cuNDA_device );

//template EXPORTGPUCORE void sqrt<float>( cuNDArray<float>*, cuNDA_device );

//template EXPORTGPUCORE void reciprocal_sqrt<float>( cuNDArray<float>*, cuNDA_device );



//template EXPORTGPUCORE void abs<float>( cuNDArray<float>*, cuNDA_device );
/*
template EXPORTGPUCORE void abs<floatd1>( cuNDArray<floatd1>*, cuNDA_device );
template EXPORTGPUCORE void abs<floatd2>( cuNDArray<floatd2>*, cuNDA_device );
template EXPORTGPUCORE void abs<floatd3>( cuNDArray<floatd3>*, cuNDA_device );
template EXPORTGPUCORE void abs<floatd4>( cuNDArray<floatd4>*, cuNDA_device );
*/
template EXPORTGPUCORE void threshold_min<float>(float, cuNDArray<float>*, float, cuNDA_device );
template EXPORTGPUCORE void threshold_min<float>(cuNDArray<float>*, cuNDArray<float>*, float, cuNDA_device );
template EXPORTGPUCORE void threshold_amin<float>(cuNDArray<float>*, cuNDArray<float>*, float, cuNDA_device );

template EXPORTGPUCORE void threshold_min<float_complext>(float_complext, cuNDArray<float_complext>*, float_complext, cuNDA_device );
template EXPORTGPUCORE void threshold_max<float>(float, cuNDArray<float>*, float, cuNDA_device );
template EXPORTGPUCORE void threshold_max<float_complext>(float_complext, cuNDArray<float_complext>*, float_complext, cuNDA_device );
template EXPORTGPUCORE void threshold_max<float>(cuNDArray<float>*, cuNDArray<float>*, float, cuNDA_device );


template EXPORTGPUCORE void rss_normalize<float>( cuNDArray<float>*, unsigned int, cuNDA_device );
template EXPORTGPUCORE void rss_normalize<float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device );


//template EXPORTGPUCORE void add<float_complext>( float_complext, cuNDArray<float_complext>*, cuNDA_device );
//template EXPORTGPUCORE void add<float>( float, cuNDArray<float>*, cuNDA_device );

//template EXPORTGPUCORE void scal<float>( float, cuNDArray<float_complext>*, cuNDA_device );

//template EXPORTGPUCORE void scale<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
//template EXPORTGPUCORE void scale<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );
//template<> EXPORTGPUCORE void scale_conj<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

//template EXPORTGPUCORE void scale<float>( cuNDArray<float>*, cuNDArray<float_complext>*, cuNDA_device );

//template EXPORTGPUCORE void axpy<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
//template EXPORTGPUCORE void axpy<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<float,1>(uintd1, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<float_complext,1>(uintd1, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<float,2>(uintd2, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<float_complext,2>(uintd2, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<float,3>(uintd3, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<float_complext,3>(uintd3, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<float,4>(uintd4, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<float_complext,4>(uintd4, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<float,float,2>(float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<float,float_complext,2>(float, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<float,float,3>(float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<float,float_complext,3>(float, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<float,float,4>(float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<float,float_complext,4>(float, cuNDArray<float_complext>*, cuNDA_device );

//template EXPORTGPUCORE float dot<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
//template EXPORTGPUCORE float_complext dot<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );
/*
template EXPORTGPUCORE float asum<float,float>( cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE float asum<float,float_complext>( cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE float nrm2<float,float>( cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE float nrm2<float,float_complext>( cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void axpy<float>( float, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void axpy<float_complext>( float_complext, cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void scal<float>( float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void scal<float_complext>( float_complext, cuNDArray<float_complext>*, cuNDA_device );
*/


template EXPORTGPUCORE void shrink1<float,float>( float, cuNDArray<float>*, cuNDArray<float>* );
template EXPORTGPUCORE void shrink1<float,float_complext>( float, cuNDArray<float_complext>*, cuNDArray<float_complext>* );

template EXPORTGPUCORE void shrinkd<float,float>( float, cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>* );
template EXPORTGPUCORE void shrinkd<float,float_complext>( float, cuNDArray<float>*, cuNDArray<float_complext>*, cuNDArray<float_complext>* );

template EXPORTGPUCORE
void origin_mirror<float,1>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
void origin_mirror<float,2>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
void origin_mirror<float,3>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
void origin_mirror<float,4>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);

template EXPORTGPUCORE void
origin_mirror<float_complext,1>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
origin_mirror<float_complext,2>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
origin_mirror<float_complext,3>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
origin_mirror<float_complext,4>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);


// Instanciation -- double precision

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
sum<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<complext< double> > >
sum<complext<double> >( cuNDArray<complext< double> >*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<1>::Type> > 
sum<doubled<1>::Type>( cuNDArray<doubled<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<2>::Type> > 
sum<doubled<2>::Type>( cuNDArray<doubled<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<3>::Type> >
 sum<doubled<3>::Type>( cuNDArray<doubled<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<4>::Type> >
 sum<doubled<4>::Type>( cuNDArray<doubled<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );



template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
expand<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
expand<double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device);


template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
abs<double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
abs<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
cNorm<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
cNorm<double,double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );


template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
minimum<double>( cuNDArray<double>*, cuNDArray<double>*,device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
maximum<double>( cuNDArray<double>*, cuNDArray<double>*,device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm<double,1>( cuNDArray<doubled<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm<double,2>( cuNDArray<doubled<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm<double,3>( cuNDArray<doubled<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm<double,4>( cuNDArray<doubled<4>::Type>*, cuNDA_device, cuNDA_device );
*/
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm_squared<double,double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
 norm_squared<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );
*/
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm_squared<double,1>( cuNDArray<doubled<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm_squared<double,2>( cuNDArray<doubled<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm_squared<double,3>( cuNDArray<doubled<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
norm_squared<double,4>( cuNDArray<doubled<4>::Type>*, cuNDA_device, cuNDA_device );


template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
ss<double,double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
ss<double,double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
ss<double_complext, double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
rss<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
rss<double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

//template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
//rss<double_complext, double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
reciprocal_rss<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
reciprocal_rss<double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );
*/
//template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
//reciprocal_rss<double_complext, double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
correlation<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
correlation<double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );

//template EXPORTGPUCORE void axpy<double>( cuNDArray<double>*, cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
crop<double,1>( uintd1, cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
crop<complext<double> ,1>( uintd1, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<complext<double> ,2>( uintd2, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<complext<double> ,3>( uintd3, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<complext<double> ,4>( uintd4, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );



template EXPORTGPUCORE void
crop<vector_td<double,1>,1>( uintd1, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,2>,1>( uintd1, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,3>,1>( uintd1, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,4>,1>( uintd1, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<double,1>,2>( uintd2, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,2>,2>( uintd2, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,3>,2>( uintd2, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,4>,2>( uintd2, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<double,1>,3>( uintd3, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,2>,3>( uintd3, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,3>,3>( uintd3, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,4>,3>( uintd3, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
crop<vector_td<double,1>,4>( uintd4, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,2>,4>( uintd4, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,3>,4>( uintd4, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
crop<vector_td<double,4>,4>( uintd4, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
expand_with_zero_fill<double,1>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<double_complext,1>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
expand_with_zero_fill<double,2>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<double_complext,2>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
expand_with_zero_fill<double,3>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<double_complext,3>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
expand_with_zero_fill<double,4>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
expand_with_zero_fill<double_complext,4>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
real_to_complext<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
complext_to_real<double>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
downsample<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
downsample<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
downsample<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
downsample<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_nn<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_nn<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_nn<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_nn<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_lin<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_lin<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_lin<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
upsample_lin<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );
/*
template EXPORTGPUCORE void clear<double>( cuNDArray<double>*,double, cuNDA_device );
template EXPORTGPUCORE void clear<double_complext>( cuNDArray<double_complext>*,double_complext, cuNDA_device );

template EXPORTGPUCORE void reciprocal<double>( cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void reciprocal<double_complext>( cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void sqrt<double>( cuNDArray<double>*, cuNDA_device );

template EXPORTGPUCORE void reciprocal_sqrt<double>( cuNDArray<double>*, cuNDA_device );

template EXPORTGPUCORE void abs<double>( cuNDArray<double>*, cuNDA_device );


template EXPORTGPUCORE void abs<doubled1>( cuNDArray<doubled1>*, cuNDA_device );
template EXPORTGPUCORE void abs<doubled2>( cuNDArray<doubled2>*, cuNDA_device );
template EXPORTGPUCORE void abs<doubled3>( cuNDArray<doubled3>*, cuNDA_device );
template EXPORTGPUCORE void abs<doubled4>( cuNDArray<doubled4>*, cuNDA_device );
*/
template EXPORTGPUCORE void threshold_min<double>(double, cuNDArray<double>*, double, cuNDA_device );
template EXPORTGPUCORE void threshold_min<double_complext>(double_complext, cuNDArray<double_complext>*, double_complext, cuNDA_device );
template EXPORTGPUCORE void threshold_min<double>(cuNDArray<double>*, cuNDArray<double>*, double, cuNDA_device );
template EXPORTGPUCORE void threshold_amin<double>(cuNDArray<double>*, cuNDArray<double>*, double, cuNDA_device );
template EXPORTGPUCORE void threshold_max<double>(double, cuNDArray<double>*, double, cuNDA_device );
template EXPORTGPUCORE void threshold_max<double_complext>(double_complext, cuNDArray<double_complext>*, double_complext, cuNDA_device );
template EXPORTGPUCORE void threshold_max<double>(cuNDArray<double>*, cuNDArray<double>*, double, cuNDA_device );

template EXPORTGPUCORE void rss_normalize<double>( cuNDArray<double>*, unsigned int, cuNDA_device );
template EXPORTGPUCORE void rss_normalize<double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device );

/*
template EXPORTGPUCORE void add<double_complext>( double_complext, cuNDArray<double_complext>*, cuNDA_device );
template EXPORTGPUCORE void add<double>( double, cuNDArray<double>*, cuNDA_device );

template EXPORTGPUCORE void scal<double>( double, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void scale<double>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void scale<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void scale<double>( cuNDArray<double>*, cuNDArray<double_complext>*, cuNDA_device );
*/
//template EXPORTGPUCORE void scale_conj<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

//template EXPORTGPUCORE void axpy<double>( cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
//template EXPORTGPUCORE void axpy<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<double,1>(uintd1, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<double_complext,1>(uintd1, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<double,2>(uintd2, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<double_complext,2>(uintd2, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<double,3>(uintd3, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<double_complext,3>(uintd3, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<double,4>(uintd4, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<double_complext,4>(uintd4, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<double,double,2>(double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<double,double_complext,2>(double, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<double,double,3>(double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<double,double_complext,3>(double, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void zero_fill_border<double,double,4>(double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void zero_fill_border<double,double_complext,4>(double, cuNDArray<double_complext>*, cuNDA_device );

/*
template EXPORTGPUCORE double dot<double>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE double_complext dot<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE double asum<double,double>( cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE double asum<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE double nrm2<double,double>( cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE double nrm2<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void axpy<double>( double, cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void axpy<double_complext>( double_complext, cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void scal<double>( double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void scal<double_complext>( double_complext, cuNDArray<double_complext>*, cuNDA_device );
*/

template EXPORTGPUCORE void shrink1<double,double>( double, cuNDArray<double>*, cuNDArray<double>* );
template EXPORTGPUCORE void shrink1<double,double_complext>( double, cuNDArray<double_complext>*, cuNDArray<double_complext>* );

template EXPORTGPUCORE void shrinkd<double,double>( double, cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>* );
template EXPORTGPUCORE void shrinkd<double,double_complext>( double, cuNDArray<double>*, cuNDArray<double_complext>*, cuNDArray<double_complext>* );

template EXPORTGPUCORE
void origin_mirror<double,1>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
void origin_mirror<double,2>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
void origin_mirror<double,3>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
void origin_mirror<double,4>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);

template EXPORTGPUCORE void
origin_mirror<double_complext,1>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
origin_mirror<double_complext,2>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
origin_mirror<double_complext,3>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
origin_mirror<double_complext,4>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> >
squaredNorm( cuNDArray<float> *, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> >
squaredNorm( cuNDArray<float_complext> *, unsigned int ,cuNDA_device , cuNDA_device);

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> >
squaredNorm( cuNDArray<double> *, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> >
squaredNorm( cuNDArray<double_complext> *, unsigned int ,cuNDA_device , cuNDA_device);


