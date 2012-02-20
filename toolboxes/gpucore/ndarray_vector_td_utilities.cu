#include "ndarray_vector_td_utilities.h"
#include "real_utilities.h"
#include "real_utilities_device.h"
#include "check_CUDA.h"

#include <cublas_v2.h>

#include <vector>
#include <cmath>


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
static bool prepare( cuNDA_device compute_device, int *cur_device, int *old_device,
		     cuNDArray<I1> *in1,       cuNDArray<I1> **in1_int, 
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> **in2_int = 0x0, 
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> **in3_int = 0x0 ) 
{
  // Test validity of D
  if( D==0 || D>3 ){
    cerr << endl << ">>>Internal error<<< :prepare: D out of range" << endl;
    return false;
  }

  if( !cur_device || !old_device ){
    cerr << endl << ">>>Internal error<<< :prepare: device ids 0x0" << endl;
    return false;
  }

  // Test validity of input pointer
  if( !in1 || !in1_int ){
    cerr << endl << "unable to process 0x0 input";
    return false;
  }
  if( D>1 && (!in2 || !in2_int) ){
    cerr << endl << "unable to process 0x0 input";
    return false;
  }
  if( D>2 && (!in3 || !in3_int) ){
    cerr << endl << "unable to process 0x0 input";
    return false;
  }
  
  // Get current Cuda device
  if( cudaGetDevice(old_device) != cudaSuccess ) {
    cerr << endl << "unable to get device no";
    return false;
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
    cerr << endl << ">>>Internal error<<< :prepare: unknown compute mode" << endl;
    return false;
  }

  if( *cur_device != *old_device && cudaSetDevice(*cur_device) != cudaSuccess) {
    cerr << endl << "unable to set device no";
    return false;
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
  
  return true;
}  

// Restore active device and free internal memory 
//
template< unsigned int D, typename I1, typename O, typename I2, typename I3 > 
static bool restore( int old_device,
		     cuNDArray<I1> *in1, cuNDArray<I1> *in1_int, 
		     unsigned int out_idx = 0, cuNDA_device alloc_device = CUNDA_NDARRAY_DEVICE, cuNDArray<O>  *out = 0x0,
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> *in2_int = 0x0, 
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> *in3_int = 0x0 )
{
  // Test validity of D
  if( D==0 || D>3 ){
    cerr << endl << ">>>Internal error<<< :prepare: D out of range" << endl;
    return false;
  }

  // Test validity of input pointer
  if( !in1 || !in1_int ){
    cerr << endl << "unable to process 0x0 input";
    return false;
  }
  if( D>1 && (!in2 || !in2_int) ){
    cerr << endl << "unable to process 0x0 input";
    return false;
  }
  if( D>2 && (!in3 || !in3_int) ){
    cerr << endl << "unable to process 0x0 input";
    return false;
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
      cerr << endl << ">>>Internal error<<< :restore: array index out of range" << endl;
      return false;
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
    cerr << endl << ">>>Internal error<<< :restore: illegal device specified" << endl;
    return false;
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
    cerr << endl << "unable to get device no";
    return false;
  }

  // Restore old device
  if( device != old_device && cudaSetDevice(old_device) != cudaSuccess) {
    cerr << endl << "unable to restore device no";
    return false;
  }
    
  return true;
}

// Initialize static variables
//
static bool initialize_static_variables()
{
  // This function is executed only once
  if( num_devices ) return true;

  if( cudaGetDeviceCount( &num_devices ) != cudaSuccess) {
    cout << endl << "Error: no Cuda devices present.";
    num_devices = 0;
    return false;
  }

  int old_device;
  if( cudaGetDevice(&old_device) != cudaSuccess ) {
    cerr << endl << "Error: unable to get device no";
    return false;
  }

  warp_size = new int[num_devices];
  max_blockdim = new int[num_devices];
  max_griddim = new int[num_devices];
  handle = new cublasHandle_t[num_devices];

  if( !warp_size || !max_blockdim || !max_griddim || !handle ) {
    cout << endl << "Error: trivial malloc failed!" << endl ;
    return false;
  }

  for( int device=0; device<num_devices; device++ ){

    if( cudaSetDevice(device) != cudaSuccess ) {
      cerr << endl << "Error: unable to set device no";
      return false;
    }
    
    cudaDeviceProp deviceProp; 
    
    if( cudaGetDeviceProperties( &deviceProp, device ) != cudaSuccess) {
      cout << endl << "Error: unable to determine device properties." << endl ;
      return false;
    }

    warp_size[device] = deviceProp.warpSize;
    max_blockdim[device] = deviceProp.maxThreadsDim[0];
    max_griddim[device] = deviceProp.maxGridSize[0];

    if (cublasCreate(&handle[device]) != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "Error: unable to create cublas handle for device " << device << std::endl;
      return false;
    }

    cublasSetPointerMode( handle[device], CUBLAS_POINTER_MODE_HOST );

  }
  
  if( cudaSetDevice(old_device) != cudaSuccess ) {
    cerr << endl << "Error: unable to restore device no";
    return false;
  }
  
  return true;
}

// Common block/grid configuration utility
//
static bool setup_grid( unsigned int cur_device, unsigned int number_of_elements, 
			dim3 *blockDim, dim3* gridDim, unsigned int num_batches=1 )
{
  if( num_devices==0 && !initialize_static_variables() ){
    cerr << endl << "system device error"; 
    return false;
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
   
  if( gridDim->x > max_griddim[cur_device] || gridDim->y > max_griddim[cur_device] )
    return false;
  else 
    return true;
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

// Norm
//
template<class REAL, unsigned int D> __global__ void 
cuNDA_norm_kernel( typename reald<REAL,D>::Type *in, REAL *out, 
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
cuNDA_norm( cuNDArray<typename reald<REAL,D>::Type> *in,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device; 
  cuNDArray< typename reald<REAL,D>::Type > *in_int;

  // Prepare
  if( !prepare<1,typename reald<REAL,D>::Type,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_norm: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_norm: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) cuNDA_norm_kernel<REAL,D><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in_int->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,typename reald<REAL,D>::Type,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_norm: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
 
  return out;
}
// cAbs
//
template<class REAL, class T> __global__ void
cuNDA_cAbs_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    T val = in[idx];
    out[idx] = abs(val);
  }
}

// Abs
//
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<REAL> > 
cuNDA_cAbs( cuNDArray<T> *in,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Prepare 
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_cAbs: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_cAbs: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) cuNDA_cAbs_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_cAbs: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}

// cNorm
//
template<class REAL, class T> __global__ void
cuNDA_cNorm_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    T val = in[idx];
    out[idx] = norm(val);
  }
}

// cNorm
//
template<class REAL, class T>
boost::shared_ptr< cuNDArray<REAL> >
cuNDA_cNorm( cuNDArray<T> *in,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Prepare
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_abs: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_abs: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get()); 
  if( out.get() != 0x0 ) cuNDA_cNorm_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
 
  // Restore 
  if( !restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_abs: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}


// Norm squared
//
template<class REAL, unsigned int D> __global__ 
void cuNDA_norm_squared_kernel( typename reald<REAL,D>::Type *in, REAL *out, unsigned int number_of_elements )
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
cuNDA_norm_squared( cuNDArray<typename reald<REAL,D>::Type> *in,
		    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device;
  cuNDArray< typename reald<REAL,D>::Type > *in_int;

  // Prepare
  if( !prepare<1,typename reald<REAL,D>::Type,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_norm_squared: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_norm_squared: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) cuNDA_norm_squared_kernel<REAL,D><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  if( !restore<1,typename reald<REAL,D>::Type,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_norm_squared: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}


// Sum
//
template<class T> __global__ void
cuNDA_sum_kernel( T *in, T *out, 
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
cuNDA_sum( cuNDArray<T> *in, unsigned int dim,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_sum: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }
  
  // Some validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_sum: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<T> >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_sum: dimension out of range." << endl; 
    return boost::shared_ptr< cuNDArray<T> >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_sum: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }
 
  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) cuNDA_sum_kernel<T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_sum: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }

  return out;
}

// Expand
//
template<class T> __global__ void
cuNDA_expand_kernel( T *in, T *out, 
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
cuNDA_expand( cuNDArray<T> *in, unsigned int new_dim_size,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_expand: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }
   
  unsigned int number_of_elements = in->get_number_of_elements()*new_dim_size;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_expand: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }
 
  // Find element stride
  vector<unsigned int> dims = *in->get_dimensions(); 
  dims.push_back(new_dim_size);

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) cuNDA_expand_kernel<T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), number_of_elements, new_dim_size );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_expand: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }

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
cuNDA_ss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _ss<REAL,T>(idx, in, stride, number_of_batches); 
  }
}

// SS
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<REAL> > 
_cuNDA_ss( cuNDArray<T> *in, unsigned int dim, 
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_ss: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_ss: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_ss: dimension out of range." << endl; 
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_ss: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if ( out.get() != 0x0 ) cuNDA_ss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  if( !restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_ss: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
cuNDA_ss<float,float>( cuNDArray<float> *in, unsigned int dim,
		       cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_ss<float, float>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
cuNDA_ss<double,double>( cuNDArray<double> *in, unsigned int dim,
			 cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_ss<double, double>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
cuNDA_ss<float,float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
				      cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_ss<float, float_complext>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
cuNDA_ss<double,double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
					cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_ss<double, double_complext>(in, dim, alloc_device, compute_device);
}

// cSS
template<class REAL, class T> __global__ void
cuNDA_css_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
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
cuNDA_css( cuNDArray<complext<REAL> > *in, unsigned int dim,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_css: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_css: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_css: dimension out of range." << endl; 
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_css: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<complext<REAL> >( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(&dims);
  if ( out.get() != 0x0 ) cuNDA_css_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  if( !restore<1,complext<REAL>,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_css: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_ss<float_complext, float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
						      cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return cuNDA_css<float>( in, dim, alloc_device, compute_device );
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_ss<double_complext, double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
							cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return cuNDA_css<double>( in, dim, alloc_device, compute_device );
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
cuNDA_rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _rss<REAL,T>(idx, in, stride, number_of_batches); 
  }
}

// RSS
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<REAL> > 
_cuNDA_rss( cuNDArray<T> *in, unsigned int dim,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_rss: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_rss: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_rss: dimension out of range." << endl; 
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_rss: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims); 
  if ( out.get() != 0x0 ) cuNDA_rss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  if( !restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_rss: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
cuNDA_rss<float,float>( cuNDArray<float> *in, unsigned int dim,
			cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_rss<float, float>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
cuNDA_rss<double,double>( cuNDArray<double> *in, unsigned int dim,
			  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_rss<double, double>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
cuNDA_rss<float,float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
				       cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_rss<float, float_complext>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
cuNDA_rss<double,double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
					 cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_rss<double, double_complext>(in, dim, alloc_device, compute_device);
}

// cRSS
template<class REAL, class T> __global__ void
cuNDA_crss_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
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
cuNDA_crss( cuNDArray<complext<REAL> > *in, unsigned int dim,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_crss: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_crss: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_crss: dimension out of range." << endl; 
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_crss: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<complext<REAL> >( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(&dims);
  if ( out.get() != 0x0 ) cuNDA_crss_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  if( !restore<1,complext<REAL>,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_crss: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
  
  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_rss<float_complext, float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
						       cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return cuNDA_crss<float>( in, dim, alloc_device, compute_device );
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_rss<double_complext, double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
							 cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return cuNDA_crss<double>( in, dim, alloc_device, compute_device );
}

// reciprocal RSS
template<class REAL, class T> __global__ void
cuNDA_reciprocal_rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = 1/(_rss<REAL,T>(idx, in, stride, number_of_batches));
  }
}

// Reciprocal RSS
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<REAL> > 
_cuNDA_reciprocal_rss( cuNDArray<T> *in, unsigned int dim, 	   
		       cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_reciprocal_rss: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_reciprocal_rss: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_reciprocal_rss: dimension out of range." << endl; 
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_reciprocal_rss: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if ( out.get() != 0x0 ) cuNDA_reciprocal_rss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  if( !restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_reciprocal_rss: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
cuNDA_reciprocal_rss<float,float>( cuNDArray<float> *in, unsigned int dim,
				   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_reciprocal_rss<float, float>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
cuNDA_reciprocal_rss<double,double>( cuNDArray<double> *in, unsigned int dim,
				     cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_reciprocal_rss<double, double>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
cuNDA_reciprocal_rss<float,float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
						  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_reciprocal_rss<float, float_complext>(in, dim, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> > 
cuNDA_reciprocal_rss<double,double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
						    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_reciprocal_rss<double, double_complext>(in, dim, alloc_device, compute_device);
}

// cReciprocal RSS
template<class REAL, class T> __global__ void
cuNDA_creciprocal_rss_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
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
cuNDA_creciprocal_rss( cuNDArray<complext<REAL> > *in, unsigned int dim,
		       cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_creciprocal_rss: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_creciprocal_rss: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_creciprocal_rss: dimension out of range." << endl; 
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_creciprocal_rss: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<complext<REAL> >( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(&dims);
  if ( out.get() != 0x0 ) cuNDA_creciprocal_rss_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  if( !restore<1,complext<REAL>,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_creciprocal_rss: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_reciprocal_rss<float_complext, float_complext>( cuNDArray<float_complext> *in, unsigned int dim,
								  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return cuNDA_creciprocal_rss<float>( in, dim, alloc_device, compute_device );
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_reciprocal_rss<double_complext, double_complext>( cuNDArray<double_complext> *in, unsigned int dim,
								    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return cuNDA_creciprocal_rss<double>( in, dim, alloc_device, compute_device );
}

// Build correlation matrix
template<class REAL, class T> __global__ void
cuNDA_correlation_kernel( T *in, T *corrm, unsigned int num_batches, unsigned int num_elements )
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
_cuNDA_correlation( cuNDArray<T> *in,
		    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_correlation: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_correlation: underdimensioned." << endl; 
    return boost::shared_ptr< cuNDArray<T> >();
  }
 
  unsigned int number_of_batches = in->get_size(in->get_number_of_dimensions()-1);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  dim3 blockDim(((max_blockdim[old_device]/number_of_batches)/warp_size[old_device])*warp_size[old_device], number_of_batches);

  if( blockDim.x == 0 ){
    cout << endl << "cuNDA_correlation: correlation dimension exceeds device capacity." << endl; 
    return boost::shared_ptr< cuNDArray<T> >();
  }
  
  dim3 gridDim((number_of_elements+blockDim.x-1)/blockDim.x);

  // Invoke kernel
  vector<unsigned int> dims = *in->get_dimensions(); dims.push_back(number_of_batches);
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) cuNDA_correlation_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_correlation: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }

  return out;
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> > 
cuNDA_correlation<float>( cuNDArray<float> *data,
			  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_correlation<float,float>(data, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> >
cuNDA_correlation<double>( cuNDArray<double> *data,
			   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_correlation<double,double>(data, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_correlation<float_complext>( cuNDArray<float_complext> *data,
					 cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_correlation<float,float_complext>(data, alloc_device, compute_device);
}

template<> EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_correlation<double_complext>( cuNDArray<double_complext> *data,
					  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  return _cuNDA_correlation<double,double_complext>(data, alloc_device, compute_device);
}

// Real to complext
template<class REAL> __global__ void
cuNDA_real_to_complext_kernel( REAL *in, complext<REAL> *out, unsigned int num_elements )
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
cuNDA_real_to_complext( cuNDArray<REAL> *in, 	    
			cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_real_to_complext: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }
 
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_real_to_complext: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) cuNDA_real_to_complext_kernel<REAL><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements());
  
  CHECK_FOR_CUDA_ERROR();
  
  // Restore
  if( !restore<1,REAL,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_real_to_complext: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  return out;
}

// complext to real by cropping the imaginary component 
template<class REAL> __global__ void
cuNDA_complext_to_real_kernel( complext<REAL> *in, REAL *out, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < num_elements ){
    out[idx] = in[idx].vec[0];
  }
}

// Convert complext to real by cropping the imaginary component 
template<class REAL>  
boost::shared_ptr< cuNDArray<REAL> > 
cuNDA_complext_to_real( cuNDArray<complext<REAL> > *in,
			cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_complext_to_real: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
 
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_real_to_complext: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());  
  if( out.get() != 0x0 ) cuNDA_complext_to_real_kernel<REAL><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements());
  
  CHECK_FOR_CUDA_ERROR();
  
  // Restore
  if( !restore<1,complext<REAL>,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_real_to_complext: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}

// Downsample
template<class REAL, unsigned int D> __global__ void
cuNDA_downsample_kernel( REAL *in, REAL *out, 
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

    REAL res = REAL(0);

    for( unsigned int i=0; i<num_adds; i++ ){
      const typename uintd<D>::Type local_co = idx_to_co<D>( i, twos );
      const unsigned int in_idx = co_to_idx<D>(co_in+local_co, matrix_size_in)+frame_offset*prod(matrix_size_in); 
      res += in[in_idx];
    }
    
    out[idx] = res/(REAL)num_adds;
  }
}

// Downsample
template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<REAL> > 
cuNDA_downsample( cuNDArray<REAL> *in,
		  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_downsample: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    cerr << endl << "cuNDA_downsample: the number of array dimensions should be at least D" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  for( unsigned int d=0; d<D; d++ ){
    if( (in->get_size(d)%2) == 1 ){
      cerr << endl << "cuNDA_downsample: uneven array dimensions not accepted" << endl;
      return boost::shared_ptr< cuNDArray<REAL> >();
    }
  }
  
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = matrix_size_in >> 1;

  unsigned int number_of_elements = prod(matrix_size_out);
  unsigned int number_of_batches = 1;

  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;

  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches ) ){
    cerr << endl << "cuNDA_scale: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out); 
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    cuNDA_downsample_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_downsample: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}

// Nearest neighbor upsampling
template<class REAL, unsigned int D> __global__ void
cuNDA_upsample_nn_kernel( REAL *in, REAL *out, 
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
cuNDA_upsample_nn( cuNDArray<REAL> *in,
		   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_upsample: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    cerr << endl << "cuNDA_upsample: the number of array dimensions should be at least D" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
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

  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches ) ){
    cerr << endl << "cuNDA_scale: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out); 
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    cuNDA_upsample_nn_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_upsample: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}

// Utility to check if all neighbors required for the linear interpolation exists
// 

template<class REAL, unsigned int D> __device__ 
bool is_border_pixel( typename uintd<D>::Type co, typename uintd<D>::Type dims )
{
  for( unsigned int dim=0; dim<D; dim++ ){
    if( co.vec[dim] == 0 || co.vec[dim] == (dims.vec[dim]-1) )
      return true;
  }
  return false;
}

// Linear upsampling
template<class REAL, unsigned int D> __global__ void
cuNDA_upsample_lin_kernel( REAL *in, REAL *out, 
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

	const typename uintd<D>::Type ones = to_vector_td<unsigned int,D>(1);
	const typename uintd<D>::Type twos = to_vector_td<unsigned int,D>(2);
	const typename uintd<D>::Type stride = idx_to_co<D>( i, twos );

	typename uintd<D>::Type co_in = ((co_out-ones)>>1)+stride;
	
	// Read corresponding pixel value
	//
	
	const unsigned int in_idx = co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in); 
	REAL value = in[in_idx];
	
	// Determine weight
	//
	
	REAL weight = REAL(1);
	
	for( unsigned int dim=0; dim<D; dim++ ){
	  
	  if( stride.vec[dim] == (co_out.vec[dim]%2) ) {
	    weight *= REAL(0.25);
	  }
	  else{
	    weight *= REAL(0.75);
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
cuNDA_upsample_lin( cuNDArray<REAL> *in,
		    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ){
    cerr << endl << "cuNDA_upsample: unable to prepare device(s)" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    cerr << endl << "cuNDA_upsample: the number of array dimensions should be at least D" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
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

  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches ) ){
    cerr << endl << "cuNDA_scale: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out); 
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    cuNDA_upsample_lin_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() ) ){
    cerr << endl << "cuNDA_upsample: unable to restore device" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  return out;
}

// Clear
template<class T> __global__ 
void cuNDA_clear_kernel( T *in_out, T val, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = val;
  }
}

// Clear
template<class T> 
bool cuNDA_clear( cuNDArray<T> *in_out, T val,
		  cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_clear: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_clear: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_clear_kernel<T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), val, in_out->get_number_of_elements() );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_clear: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Abs
template<class T> __global__ 
void cuNDA_abs_kernel( T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T val = in_out[idx]; 
    in_out[idx] = abs(val);
  }
}

// Abs
template<class T>  
bool cuNDA_abs( cuNDArray<T> *in_out,
		cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_abs: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_abs: block/grid configuration out of range" << endl;
    return false;
  }
 
  // Invoke kernel
  cuNDA_abs_kernel<T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), in_out->get_number_of_elements() );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_abs: unable to restore device" << endl;
    return false;
  }

  return true;
}
// threshold
template<class T> __global__
void cuNDA_threshold_min_kernel(T min,  T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    if (in_out[idx] < min) in_out[idx] = T(0);
  }
}


//Threshold
template<class T>
bool cuNDA_threshold_min(T min, cuNDArray<T> *in_out,
		cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_threshold_min: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_threshold_min: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_threshold_min_kernel<T><<< gridDim, blockDim >>>(min, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_threshold_min: unable to restore device" << endl;
    return false;
  }

  return true;
}

// threshold
template<class T> __global__
void cuNDA_threshold_max_kernel(T max,  T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    if (in_out[idx] > max) in_out[idx] = T(0);
  }
}


//Threshold
template<class T>
bool cuNDA_threshold_max(T max, cuNDArray<T> *in_out,
		cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_threshold_min: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_threshold_min: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_threshold_max_kernel<T><<< gridDim, blockDim >>>(max, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_threshold_min: unable to restore device" << endl;
    return false;
  }

  return true;
}


// Reciprocal
template<class T> __global__ 
void cuNDA_reciprocal_kernel( T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = 1/(in_out[idx]);
  }
}

// Reciprocal
template<class T> 
bool cuNDA_reciprocal( cuNDArray<T> *in_out,
		       cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_reciprocal: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_reciprocal: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_reciprocal_kernel<T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), in_out->get_number_of_elements() );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_reciprocal: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Square root
template<class T> __global__ 
void cuNDA_sqrt_kernel( T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = sqrt(in_out[idx]);
  }
}

// Square root
template<class T> 
bool cuNDA_sqrt( cuNDArray<T> *in_out,
		 cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_sqrt: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_sqrt: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_sqrt_kernel<T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), in_out->get_number_of_elements() );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_sqrt: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Reciprocal square root
template<class T> __global__ 
void cuNDA_reciprocal_sqrt_kernel( T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = gad_rsqrt(in_out[idx]);
  }
}

// Square root
template<class T> 
bool cuNDA_reciprocal_sqrt( cuNDArray<T> *in_out,
			    cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_reciprocal_sqrt: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_reciprocal_sqrt: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_reciprocal_sqrt_kernel<T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), in_out->get_number_of_elements() );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_reciprocal_sqrt: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Normalized RSS
template<class REAL, class T> __global__ void
cuNDA_rss_normalize_kernel( T *in_out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
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
template<class REAL, class T> 
bool cuNDA_rss_normalize( cuNDArray<T> *in_out, unsigned int dim,
			  cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_rss_normalize: unable to prepare device(s)" << endl;
    return false;
  }

  // Validity checks
  if( !(in_out->get_number_of_dimensions()>1) ){
    cout << endl << "cuNDA_rss_normalize: underdimensioned." << endl; 
    return false;
  }
 
  if( dim > in_out->get_number_of_dimensions()-1 ){
    cout << endl << "cuNDA_rss_normalize: dimension out of range." << endl; 
    return false;
  }

  unsigned int number_of_batches = in_out->get_size(dim);
  unsigned int number_of_elements = in_out->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_rss_normalize: block/grid configuration out of range" << endl;
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Find element stride
  unsigned int stride; vector<unsigned int> dims;
  find_stride<T>( in_out, dim, &stride, &dims );

  // Invoke kernel
  cuNDA_rss_normalize_kernel<REAL,T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_rss_normalize: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Add

template<class T> __global__ 
void cuNDA_add_kernel( T a, T *x, unsigned int number_of_elements )

{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    x[idx] += a;
  }
}

// Add

template<class T> 
bool cuNDA_add( T a, cuNDArray<T> *in_out,
		  cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;

  cuNDArray<T > *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_add: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_add: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel

  cuNDA_add_kernel<T><<< gridDim, blockDim >>> ( a, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );

 
  CHECK_FOR_CUDA_ERROR();

  // Restore

  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){

    cerr << endl << "cuNDA_add: unable to restore device" << endl;
    return false;
  }

  return true;
}




// Scale
template<class REAL> __global__ 
void cuNDA_scale1_kernel( REAL a, complext<REAL> *x, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    complext<REAL> in = x[idx];
    in = a*in;
    x[idx] = in;
  }
}

// Scale 
template<class REAL> 
bool cuNDA_scale( REAL a, cuNDArray<complext<REAL> > *in_out,
		  cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_scale: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_scale: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_scale1_kernel<REAL><<< gridDim, blockDim >>> ( a, in_out_int->get_data_ptr(), in_out->get_number_of_elements() );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,complext<REAL>,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_scale: unable to restore device" << endl;
    return false;
  }

  return true;
}




// Scale
template<class S, class T> __global__ 
void cuNDA_scale2_kernel( S *a, T *x, unsigned int number_of_elements, unsigned int number_of_batches )
{
   const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements*number_of_batches ){
    unsigned int frame_offset = idx/number_of_elements;
    S in_a = a[idx-frame_offset*number_of_elements];
    T in_x = x[idx];
    x[idx] = in_a*in_x;
  }
}

// Scale 
template<class T> 
bool cuNDA_scale( cuNDArray<T> *a, cuNDArray<T> *x,
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
  if( !prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, a, &a_int, x, &x_int ) ){
    cerr << endl << "cuNDA_scale: unable to prepare device(s)" << endl;
    return false;
  }

  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, num_batches ) ){
    cerr << endl << "cuNDA_scale: block/grid configuration out of range" << endl;
    return false;
  }
  
  // Invoke kernel
  cuNDA_scale2_kernel<T,T><<< gridDim, blockDim >>> ( a_int->get_data_ptr(), x_int->get_data_ptr(), number_of_elements, num_batches );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<2,T,dummy,T,dummy>( old_device, a, a_int, 2, compute_device, 0x0, x, x_int ) ){
    cerr << endl << "cuNDA_scale: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Scale 
template<class REAL> 
bool cuNDA_scale( cuNDArray<REAL> *a, cuNDArray<complext<REAL> > *x,
		  cuNDA_device compute_device )
{
  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "image dimensions mismatch, cannot scale" << endl;
    return false;
  }
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *a_int; cuNDArray<complext<REAL> > *x_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<2,REAL,complext<REAL>,dummy>( compute_device, &cur_device, &old_device, a, &a_int, x, &x_int ) ){
    cerr << endl << "cuNDA_scale: unable to prepare device(s)" << endl;
    return false;
  }

  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, num_batches ) ){
    cerr << endl << "cuNDA_scale: block/grid configuration out of range" << endl;
    return false;
  }
 
  // Invoke kernel
  cuNDA_scale2_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>> ( a_int->get_data_ptr(), x_int->get_data_ptr(), number_of_elements, num_batches );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<2,REAL, dummy, complext<REAL>,dummy>( old_device, a, a_int, 2, compute_device, 0x0, x, x_int ) ){
    cerr << endl << "cuNDA_scale: unable to restore device" << endl;
    return false;
  }

  return true;
}

// 'axpy'
template<class S, class T> __global__ 
void cuNDA_axpy_kernel( S *a, T *x, T *y, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    S in_a = a[idx];
    for( unsigned int batch=0; batch<number_of_batches; batch++ ){
      unsigned int iidx = batch*number_of_elements + idx;
      T in_x = x[iidx];
      T in_y = y[iidx];
      in_y += in_a*in_x;
      y[iidx] = in_y;
    }
  }
}



// Scale conjugate w. non conjugate
template<class S, class T> __global__ 
void cuNDA_scale_conj_kernel( S *a, T *x, unsigned int number_of_elements, unsigned int number_of_batches )
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
bool cuNDA_scale_conj( cuNDArray<T> *a, cuNDArray<T> *x,
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
  if( !prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, a, &a_int, x, &x_int ) ){
    cerr << endl << "cuNDA_scale: unable to prepare device(s)" << endl;
    return false;
  }

  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, num_batches ) ){
    cerr << endl << "cuNDA_scale: block/grid configuration out of range" << endl;
    return false;
  }
  
  // Invoke kernel
  cuNDA_scale_conj_kernel<T,T><<< gridDim, blockDim >>> ( a_int->get_data_ptr(), x_int->get_data_ptr(), number_of_elements, num_batches );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<2,T,dummy,T,dummy>( old_device, a, a_int, 2, compute_device, 0x0, x, x_int ) ){
    cerr << endl << "cuNDA_scale: unable to restore device" << endl;
    return false;
  }

  return true;
}




// '.axpy' 
template<class T> 
bool cuNDA_axpy( cuNDArray<T> *a, cuNDArray<T> *x, cuNDArray<T> *y,
		 cuNDA_device compute_device )
{
  if( x->get_number_of_elements() != y->get_number_of_elements() ){
    cout << endl << "cuNDA_axpy: image dimensions mismatch in axpy" << endl;
    return false;
  }

  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "cuNDA_axpy: image dimensions mismatch" << endl;
    return false;
  }
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *a_int, *x_int, *y_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<3,T,T,T>( compute_device, &cur_device, &old_device, a, &a_int, x, &x_int, y, &y_int ) ){
    cerr << endl << "cuNDA_axpy: unable to prepare device(s)" << endl;
    return false;
  }

  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = x->get_number_of_elements() / a->get_number_of_elements();

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_axpy: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_axpy_kernel<T,T><<< gridDim, blockDim >>> ( a_int->get_data_ptr(), x_int->get_data_ptr(), y_int->get_data_ptr(), num_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<3,T,dummy,T,T>( old_device, a, a_int, 3, compute_device, 0x0, x, x_int, y, y_int ) ){
    cerr << endl << "cuNDA_axpy: unable to restore device" << endl;
    return false;
  }

  return true;
}

// '.axpy' 
template<class REAL> 
bool cuNDA_axpy( cuNDArray<REAL> *a, cuNDArray<complext<REAL> > *x, cuNDArray<complext<REAL> > *y,
		 cuNDA_device compute_device )
{
  if( x->get_number_of_elements() != y->get_number_of_elements() ){
    cout << endl << "cuNDA_axpy: image dimensions mismatch" << endl;
    return false;
  }

  if( x->get_number_of_elements() < a->get_number_of_elements() ||
      x->get_number_of_elements() % a->get_number_of_elements() ){
    cout << endl << "cuNDA_axpy: image dimensions mismatch" << endl;
    return false;
  }
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *a_int; cuNDArray<complext<REAL> > *x_int, *y_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<3,REAL,complext<REAL>,complext<REAL> >( compute_device, &cur_device, &old_device, a, &a_int, x, &x_int, y, &y_int ) ){
    cerr << endl << "cuNDA_axpy: unable to prepare device(s)" << endl;
    return false;
  }

  unsigned int number_of_elements = a->get_number_of_elements();
  unsigned int num_batches = y->get_number_of_elements() / a->get_number_of_elements();

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, number_of_elements, &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_axpy: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_axpy_kernel<REAL, complext<REAL> ><<< gridDim, blockDim >>> ( a_int->get_data_ptr(), x_int->get_data_ptr(), y_int->get_data_ptr(), num_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<3,REAL,dummy,complext<REAL>,complext<REAL> >( old_device, a, a_int, 3, compute_device, 0x0, x, x_int, y, y_int ) ){
    cerr << endl << "cuNDA_scale: unable to restore device" << endl;
    return false;
  }

  return true;
}

template<class T> T _dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, int device );

template<> float_complext
_dot<float_complext>( cuNDArray<float_complext>* arr1, cuNDArray<float_complext>* arr2, int device )
{
  float_complext ret;
  
  if (cublasCdotc( handle[device], arr1->get_number_of_elements(),
		   (const cuComplex*) arr1->get_data_ptr(), 1, 
		   (const cuComplex*) arr2->get_data_ptr(), 1,
		   (cuComplex*) &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
    }
  
  cudaThreadSynchronize();  
  return ret;
}

template<> float 
_dot<float>( cuNDArray<float>* arr1, cuNDArray<float>* arr2, int device )
{
  float ret;

  if( cublasSdot(handle[device], arr1->get_number_of_elements(),
		 arr1->get_data_ptr(), 1, 
		 arr2->get_data_ptr(), 1,
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
    }
  
  cudaThreadSynchronize();
  return ret;
}

template<> double_complext
_dot<double_complext>( cuNDArray<double_complext>* arr1, cuNDArray<double_complext>* arr2, int device )
{
  double_complext ret;
  
  if( cublasZdotc(handle[device], arr1->get_number_of_elements(),
		  (const cuDoubleComplex*) arr1->get_data_ptr(), 1, 
		  (const cuDoubleComplex*) arr2->get_data_ptr(), 1,
		  (cuDoubleComplex*) &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
    }
  
  cudaThreadSynchronize();
  return ret;
}

template<> double 
_dot<double>( cuNDArray<double>* arr1, cuNDArray<double>* arr2, int device )
{
  double ret;

  if( cublasDdot(handle[device], arr1->get_number_of_elements(),
		 arr1->get_data_ptr(), 1, 
		 arr2->get_data_ptr(), 1,
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_dot: inner product calculation using cublas failed" << std::endl;
    }
  
  cudaThreadSynchronize();
  return ret;
}

template<class T> T
cuNDA_dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, cuNDA_device compute_device )
{
  if (arr1->get_number_of_elements() != arr2->get_number_of_elements()) {
    cout << "cuNDA_dot: array dimensions mismatch" << std::endl;
    return T(0);
  }

  if( !initialize_static_variables() ){
    cout << "cuNDA_dot: initialization failed" << std::endl;
    return T(0);
  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *arr1_int, *arr2_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, arr1, &arr1_int, arr2, &arr2_int ) ){
    cerr << endl << "cuNDA_dot: unable to prepare device(s)" << endl;
    return T(0);
  }

  T ret = _dot<T>( arr1_int, arr2_int, cur_device );  

  // Restore
  if( !restore<2,T,dummy,T,dummy>( old_device, arr1, arr1_int, 0, compute_device, 0x0, arr2, arr2_int ) ){
    cerr << endl << "cuNDA_dot: unable to restore device" << endl;
    return T(0);
  }

  return ret;
}

template<class REAL, class T> REAL _sum( cuNDArray<T>* arr, int device );

template<> float
_sum<float,float_complext>( cuNDArray<float_complext>* arr, int device )
{
  float ret;
  
  if (cublasScasum( handle[device], arr->get_number_of_elements(),
		   (const cuComplex*) arr->get_data_ptr(), 1,
		   &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
      return 0;
    }
  
  cudaThreadSynchronize();
  return ret;
}

template<> float 
_sum<float,float>( cuNDArray<float>* arr, int device )
{
  float ret;
  if( cublasSasum(handle[device], arr->get_number_of_elements(),
		 arr->get_data_ptr(), 1, 
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
      return 0;
    }
  
  cudaThreadSynchronize();
  return ret;
}

template<> double
_sum<double,double_complext>( cuNDArray<double_complext>* arr, int device )
{
  double ret;
  
  if( cublasDzasum(handle[device], arr->get_number_of_elements(),
		  (const cuDoubleComplex*) arr->get_data_ptr(), 1, 
		  &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
      return 0;
    }
  
  cudaThreadSynchronize();
  return ret;
}

template<> double 
_sum<double,double>( cuNDArray<double>* arr, int device )
{
  double ret;
  if( cublasDasum(handle[device], arr->get_number_of_elements(),
		 arr->get_data_ptr(), 1, 
		 &ret) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_sum: sum calculation using cublas failed" << std::endl;
      return 0;
    }
  
  cudaThreadSynchronize();
  return ret;
}

template<class REAL, class T> REAL
cuNDA_asum( cuNDArray<T>* arr, cuNDA_device compute_device )
{
  if( !initialize_static_variables() ){
    cout << "cuNDA_asum: initialization failed" << std::endl;
    return REAL(0);
  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *arr_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, arr, &arr_int ) ){
    cerr << endl << "cuNDA_asum: unable to prepare device(s)" << endl;
    return REAL(0);
  }

  REAL ret = _sum<REAL,T>( arr_int, cur_device );  

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, arr, arr_int, 0, compute_device ) ){
    cerr << endl << "cuNDA_asum: unable to restore device" << endl;
    return REAL(0);
  }

  return ret;
}

template<class T> bool 
_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y, int device );

template<> bool 
_axpy( float_complext a, cuNDArray<float_complext>* x, cuNDArray<float_complext>* y, int device )
{
  if( cublasCaxpy( handle[device], x->get_number_of_elements(), (cuComplex*) &a,
		   (cuComplex*) x->get_data_ptr(), 1, 
		   (cuComplex*) y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_axpy: axpy calculation using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<> bool
_axpy( float a, cuNDArray<float>* x, cuNDArray<float>* y, int device )
{
  if( cublasSaxpy( handle[device], x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1, 
		   y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_axpy: axpy calculation using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<> bool 
_axpy( double_complext a, cuNDArray<double_complext>* x, cuNDArray<double_complext>* y, int device )
{
  if( cublasZaxpy( handle[device], x->get_number_of_elements(), (cuDoubleComplex*) &a,
		   (cuDoubleComplex*) x->get_data_ptr(), 1, 
		   (cuDoubleComplex*) y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_axpy: axpy calculation using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<> bool
_axpy( double a, cuNDArray<double>* x, cuNDArray<double>* y, int device )
{
  if( cublasDaxpy( handle[device], x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1, 
		   y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_axpy: axpy calculation using cublas failed" << std::endl;
      return false;
    }
  
  return true;
} 

template<class T> bool
cuNDA_axpy( T a, cuNDArray<T>* x, cuNDArray<T>* y, cuNDA_device compute_device )
{
  if ( !x || !y || x->get_number_of_elements() != y->get_number_of_elements()) {
    cout << "cuNDA_axpy: axpy array dimensions mismatch" << std::endl;
    return false;
  }
  
   if( !initialize_static_variables() ){
    cout << "cuNDA_axpy: initialization failed" << std::endl;
    return false;
  }

 // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *x_int, *y_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, x, &x_int, y, &y_int ) ){
    cerr << endl << "cuNDA_axpy: unable to prepare device(s)" << endl;
    return false;
  }

  bool ret = _axpy<T>( a, x_int, y_int, cur_device );

  // Restore
  if( !restore<2,T,dummy,T,dummy>( old_device, x, x_int, 2, compute_device, 0x0, y, y_int ) ){
    cerr << endl << "cuNDA_axpy: unable to restore device" << endl;
    return false;
  }
  
  return ret;
}

template<class T> bool 
_scal( T a, cuNDArray<T>* x, int device );

template<> bool
_scal( float_complext a, cuNDArray<float_complext>* x, int device)
{
  if( cublasCscal( handle[device], x->get_number_of_elements(), (cuComplex*) &a,
		   (cuComplex*) x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_scal: calculation using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<> bool
_scal( float a, cuNDArray<float>* x, int device ) 
{
  if( cublasSscal( handle[device], x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_scal: calculation using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<> bool
_scal( double_complext a, cuNDArray<double_complext>* x, int device)
{
  if( cublasZscal( handle[device], x->get_number_of_elements(), (cuDoubleComplex*) &a,
		   (cuDoubleComplex*) x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_scal: calculation using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<> bool
_scal( double a, cuNDArray<double>* x, int device ) 
{
  if( cublasDscal( handle[device], x->get_number_of_elements(), &a,
		   x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS ) 
    {
      cout << "cuNDA_scal: calculation using cublas failed" << std::endl;
      return false;
    }

  return true;
}

template<class T> bool
cuNDA_scal( T a, cuNDArray<T>* x, cuNDA_device compute_device )
{
  if( !initialize_static_variables() ){
    cout << "cuNDA_scal: initialization failed" << std::endl;
    return false;
  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *x_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, x, &x_int ) ){
    cerr << endl << "cuNDA_scal: unable to prepare device(s)" << endl;
    return false;
  }

  bool ret = _scal<T>( a, x_int, cur_device );

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, x, x_int, 1, compute_device ) ){
    cerr << endl << "cuNDA_scal: unable to restore device" << endl;
    return false;
  }

  return ret;
}

// Normalize (float)
template<> EXPORTGPUCORE
float cuNDA_normalize<float>( cuNDArray<float> *data, float new_max, cuNDA_device compute_device )
{
  if( !initialize_static_variables() ){
    cout << "cuNDA_normalize: initialization failed" << std::endl;
    return 0;
  }

  unsigned int number_of_elements = data->get_number_of_elements();

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<float> *data_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,float,dummy,dummy>( compute_device, &cur_device, &old_device, data, &data_int ) ){
    cerr << endl << "cuNDA_normalize: unable to prepare device(s)" << endl;
    return 0;
  }

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
  if( !restore<1,float,dummy,dummy,dummy>( old_device, data, data_int, 1, compute_device ) ){
    cerr << endl << "cuNDA_normalize: unable to restore device" << endl;
    return 0;
  }

  CHECK_FOR_CUDA_ERROR();
  return scale;
}

// Normalize (double)
template<> EXPORTGPUCORE
double cuNDA_normalize<double>( cuNDArray<double> *data, double new_max, cuNDA_device compute_device )
{
  if( !initialize_static_variables() ){
    cout << "cuNDA_normalize: initialization failed" << std::endl;
    return 0;
  }

  unsigned int number_of_elements = data->get_number_of_elements();

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<double> *data_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,double,dummy,dummy>( compute_device, &cur_device, &old_device, data, &data_int ) ){
    cerr << endl << "cuNDA_normalize: unable to prepare device(s)" << endl;
    return 0;
  }

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
  if( !restore<1,double,dummy,dummy,dummy>( old_device, data, data_int, 1, compute_device ) ){
    cerr << endl << "cuNDA_normalize: unable to restore device" << endl;
    return 0;
  }
  
  CHECK_FOR_CUDA_ERROR();
  return scale;
}

// Crop
template<class T, unsigned int D> __global__ void
cuNDA_crop_kernel( typename uintd<D>::Type offset, typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
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
template<class T, unsigned int D> 
bool cuNDA_crop( typename uintd<D>::Type offset, 
	    cuNDArray<T> *in, cuNDArray<T> *out,
	    cuNDA_device compute_device )
{
  if( in == 0x0 || out == 0x0 ){
    cout << endl << "cuNDA_crop: 0x0 ndarray provided" << endl;
    return false;
  }

  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "cuNDA_crop: image dimensions mismatch" << endl;
    return false;
  }

  if( in->get_number_of_dimensions() < D ){
    cout << endl << "cuNDA_crop: number of image dimensions should be at least " << D << endl;
    return false;
  }

  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *out->get_dimensions() );
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }

  if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
    cout << endl << "cuNDA_crop: cropping size mismatch" << endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int ) ){
    cerr << endl << "cuNDA_crop: unable to prepare device(s)" << endl;
    return false;
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim, number_of_batches ) ){
    cerr << endl << "cuNDA_crop: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_crop_kernel<T,D><<< gridDim, blockDim >>> 
    ( offset, matrix_size_in, matrix_size_out, in_int->get_data_ptr(), out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int ) ){
    cerr << endl << "cuNDA_crop: unable to restore device" << endl;
    return false;
  }

  return out;
}

// Expand and zero fill
template<class T, unsigned int D> __global__ void
cuNDA_expand_with_zero_fill_kernel( typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
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
bool cuNDA_expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out,
				  cuNDA_device compute_device )
{ 
  if( in == 0x0 || out == 0x0 ){
    cout << endl << "cuNDA_zero_fill: 0x0 ndarray provided" << endl;
    return false;
  }

  if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
    cout << endl << "cuNDA_zero_fill: image dimensions mismatch" << endl;
    return false;
  }

  if( in->get_number_of_dimensions() < D ){
    cout << endl << "cuNDA_zero_fill: number of image dimensions should be at least " << D << endl;
    return false;
  }

  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *out->get_dimensions() );
  
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }

  if( weak_greater(matrix_size_in,matrix_size_out) ){
    cout << endl << "cuNDA_expand: size mismatch, cannot expand" << endl;
    return false;
  }
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int ) ){
    cerr << endl << "cuNDA_expand: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim, number_of_batches ) ){
    cerr << endl << "cuNDA_expand: block/grid configuration out of range" << endl;
    return false;
  }
 
  // Invoke kernel
  cuNDA_expand_with_zero_fill_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, in_int->get_data_ptr(), out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int ) ){
    cerr << endl << "cuNDA_expand_with_zero_fill: unable to restore device" << endl;
    return false;
  }

  return out;
}

// Zero fill border (rectangular)
template<class T, unsigned int D> __global__ void
cuNDA_zero_fill_border_kernel( typename uintd<D>::Type matrix_size_in, typename uintd<D>::Type matrix_size_out, 
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
bool cuNDA_zero_fill_border( typename uintd<D>::Type matrix_size_in, cuNDArray<T> *in_out,
			     cuNDA_device compute_device )
{ 
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *in_out->get_dimensions() );
 
  if( weak_greater(matrix_size_in, matrix_size_out) ){
    cout << endl << "cuNDA_zero_fill: size mismatch, cannot zero fill" << endl;
    return false;
  }
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
    number_of_batches *= in_out->get_size(d);
  }

 // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_zero_fill: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_zero_fill: block/grid configuration out of range" << endl;
    return false;
  }
 
  // Invoke kernel
  cuNDA_zero_fill_border_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, in_out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_zero_fill_border: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Zero fill border (circular)
template<class REAL, class T, unsigned int D> __global__ void
cuNDA_zero_fill_border_kernel( typename reald<REAL,D>::Type radius, typename uintd<D>::Type matrix_size_out, 
			       T *image, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx < number_of_elements ){
    
    typename reald<REAL,D>::Type half_matrix_size_out_real = to_reald<REAL,unsigned int,D>( matrix_size_out>>1 );

    const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
    typename reald<REAL,D>::Type co_out_real = to_reald<REAL,unsigned int,D>( co_out );
    
    typename reald<REAL,D>::Type co_f = abs( co_out_real - half_matrix_size_out_real );
    
    if( co_f<radius )
      ; // do nothing
    else{
      T zero = T(0);
      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
	image[idx+batch*number_of_elements] = zero;
      }
    } 
  }
}

// Zero fill border (circular)
template<class REAL, class T, unsigned int D> 
bool cuNDA_zero_fill_border( typename reald<REAL,D>::Type radius, cuNDArray<T> *in_out,
			     cuNDA_device compute_device )
{
  if( in_out->get_number_of_dimensions() < D ){
    cout << endl << "Image dimensions mismatch, cannot zero fill" << endl;
    return false;
  }
 
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *in_out->get_dimensions() );
  typename reald<REAL,D>::Type matrix_size_out_real = to_reald<REAL,unsigned int,D>( matrix_size_out );

  if( weak_greater(radius, matrix_size_out_real) ){
    cout << endl << "cuNDA_zero_fillborder:: radius extends matrix size, cannot zero fill" << endl;
    return false;
  }
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
    number_of_batches *= in_out->get_size(d);
  }

 // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ){
    cerr << endl << "cuNDA_zero_fill: unable to prepare device(s)" << endl;
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_zero_fill: block/grid configuration out of range" << endl;
    return false;
  }
 
  // Invoke kernel
  cuNDA_zero_fill_border_kernel<REAL,T,D><<< gridDim, blockDim >>> ( radius, matrix_size_out, in_out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 ) ){
    cerr << endl << "cuNDA_zero_fill_border: unable to restore device" << endl;
    return false;
  }

  return true;
}

// Shrinkage
//

template<class REAL, class T> __global__ void 
cuNDA_shrink1_kernel( REAL gamma, T *in, T *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T in_val = in[idx]; 
    REAL in_norm = norm<REAL>(in_val);
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
bool cuNDA_shrink1( REAL gamma, cuNDArray<T> *in, cuNDArray<T> *out )
{
  // TODO: multi-device handling

  if( !in || !out ){
    cerr << endl << "cuNDA_shrink1: 0x0 arrays not accepted" << endl;
    return false;
  }

  if( in->get_number_of_elements() != out->get_number_of_elements() ){
    cerr << endl << "cuNDA_shrink1: i/o arrays must have an identical number of elements" << endl;
    return false;
  }
  
  // Get current Cuda device
  int cur_device;
  if( cudaGetDevice(&cur_device) != cudaSuccess ) {
    cerr << endl << "cuNDA_shrink1 : unable to get device no";
    return false;
  }


  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_shrink1: block/grid configuration out of range" << endl;
    return false;
  }
  
  // Invoke kernel
  cuNDA_shrink1_kernel<REAL,T><<< gridDim, blockDim >>>( gamma, in->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
  
  return true;
}

template<class REAL, class T> __global__ void 
cuNDA_shrinkd_kernel( REAL gamma, REAL *s_k, T *in, T *out, unsigned int number_of_elements )
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
bool cuNDA_shrinkd( REAL gamma, cuNDArray<REAL> *s_k, cuNDArray<T> *in, cuNDArray<T> *out )
{
  // TODO: multi-device handling

  if( !in || !out || !s_k ){
    cerr << endl << "cuNDA_shrinkd: 0x0 arrays not accepted" << endl;
    return false;
  }

  if( in->get_number_of_elements() != out->get_number_of_elements() ){
    cerr << endl << "cuNDA_shrinkd: i/o arrays must have an identical number of elements" << endl;
    return false;
  }

  if( in->get_number_of_elements() != s_k->get_number_of_elements() ){
    cerr << endl << "cuNDA_shrinkd: i/o arrays must have an identical number of elements" << endl;
    return false;
  }
  
  // Get current Cuda device
  int cur_device;
  if( cudaGetDevice(&cur_device) != cudaSuccess ) {
    cerr << endl << "cuNDA_shrinkd : unable to get device no";
    return false;
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_shrinkd: block/grid configuration out of range" << endl;
    return false;
  }
  
  // Invoke kernel
  cuNDA_shrinkd_kernel<REAL,T><<< gridDim, blockDim >>>( gamma, s_k->get_data_ptr(), in->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
  
  return true;
}

// Crop
template<class T, unsigned int D> __global__ void
cuNDA_origin_mirror_kernel( typename uintd<D>::Type matrix_size, typename uintd<D>::Type origin, T *in, T *out, bool zero_fill )
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
bool cuNDA_origin_mirror( cuNDArray<T> *in, cuNDArray<T> *out, bool zero_fill, cuNDA_device compute_device )
{
  if( in == 0x0 || out == 0x0 ){
    cout << endl << "cuNDA_origin_mirror: 0x0 ndarray provided" << endl;
    return false;
  }

  if( !in->dimensions_equal(out) ){
    cout << endl << "cuNDA_origin_mirror: image dimensions mismatch" << endl;
    return false;
  }
  
  if( in->get_number_of_dimensions() != D ){
    cout << endl << "cuNDA_origin_mirror: number of image dimensions is not " << D << endl;
    return false;
  }

  typename uintd<D>::Type matrix_size = vector_to_uintd<D>( *in->get_dimensions() );
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  if( !prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int ) ){
    cerr << endl << "cuNDA_origin_mirror: unable to prepare device(s)" << endl;
    return false;
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( cur_device, prod(matrix_size), &blockDim, &gridDim ) ){
    cerr << endl << "cuNDA_origin_mirror: block/grid configuration out of range" << endl;
    return false;
  }

  // Invoke kernel
  cuNDA_origin_mirror_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size, matrix_size>>1, in_int->get_data_ptr(), out_int->get_data_ptr(), zero_fill );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  if( !restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int ) ){
    cerr << endl << "cuNDA_origin_mirror: unable to restore device" << endl;
    return false;
  }

  return out;
}

//
// Instantiation
//

// A few functions have integer support

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<int> > 
cuNDA_sum<int>( cuNDArray<int>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<1>::Type> > 
cuNDA_sum<intd<1>::Type >( cuNDArray<intd<1>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<2>::Type> > 
cuNDA_sum<intd<2>::Type >( cuNDArray<intd<2>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<3>::Type> > 
cuNDA_sum<intd<3>::Type >( cuNDArray<intd<3>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<4>::Type> > 
cuNDA_sum<intd<4>::Type >( cuNDArray<intd<4>::Type >*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<unsigned int> > 
cuNDA_sum<unsigned int>( cuNDArray<unsigned int>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<1>::Type> > 
cuNDA_sum<uintd<1>::Type>( cuNDArray<uintd<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<2>::Type> > 
cuNDA_sum<uintd<2>::Type>( cuNDArray<uintd<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<3>::Type> > 
cuNDA_sum<uintd<3>::Type>( cuNDArray<uintd<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<4>::Type> > 
cuNDA_sum<uintd<4>::Type>( cuNDArray<uintd<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<int,1>( uintd1::Type, cuNDArray<int>*, cuNDArray<int>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,1>,1>( uintd1::Type, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,2>,1>( uintd1::Type, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,3>,1>( uintd1::Type, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,4>,1>( uintd1::Type, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,1>,2>( uintd2::Type, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,2>,2>( uintd2::Type, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,3>,2>( uintd2::Type, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,4>,2>( uintd2::Type, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,1>,3>( uintd3::Type, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,2>,3>( uintd3::Type, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,3>,3>( uintd3::Type, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,4>,3>( uintd3::Type, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,1>,4>( uintd4::Type, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,2>,4>( uintd4::Type, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,3>,4>( uintd4::Type, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<int,4>,4>( uintd4::Type, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<unsigned int,1>( uintd1::Type, cuNDArray<unsigned int>*, cuNDArray<unsigned int>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,1>,1>( uintd1::Type, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,2>,1>( uintd1::Type, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,3>,1>( uintd1::Type, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,4>,1>( uintd1::Type, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,1>,2>( uintd2::Type, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,2>,2>( uintd2::Type, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,3>,2>( uintd2::Type, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,4>,2>( uintd2::Type, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,1>,3>( uintd3::Type, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,2>,3>( uintd3::Type, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,3>,3>( uintd3::Type, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,4>,3>( uintd3::Type, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,1>,4>( uintd4::Type, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,2>,4>( uintd4::Type, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,3>,4>( uintd4::Type, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<unsigned int,4>,4>( uintd4::Type, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

// Instanciation -- single precision

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_sum<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<complext<float> > >
cuNDA_sum<complext<float> >( cuNDArray<complext<float> >*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<1>::Type> > 
cuNDA_sum<floatd<1>::Type>( cuNDArray<floatd<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<2>::Type> > 
cuNDA_sum<floatd<2>::Type>( cuNDArray<floatd<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<3>::Type> >
 cuNDA_sum<floatd<3>::Type>( cuNDArray<floatd<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<4>::Type> >
 cuNDA_sum<floatd<4>::Type>( cuNDArray<floatd<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_expand<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_expand<float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device);



template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_cAbs<float,float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
cuNDA_cAbs<float,float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
cuNDA_cNorm<float,float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
cuNDA_cNorm<float,float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm<float,1>( cuNDArray<floatd<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm<float,2>( cuNDArray<floatd<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm<float,3>( cuNDArray<floatd<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm<float,4>( cuNDArray<floatd<4>::Type>*, cuNDA_device, cuNDA_device );

/*template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
cuNDA_norm_squared<float,float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
 cuNDA_norm_squared<float,float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );*/

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm_squared<float,1>( cuNDArray<floatd<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm_squared<float,2>( cuNDArray<floatd<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm_squared<float,3>( cuNDArray<floatd<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_norm_squared<float,4>( cuNDArray<floatd<4>::Type>*, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_ss<float,float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_ss<float,float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_ss<float_complext, float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_rss<float,float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_rss<float,float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_rss<float_complext, float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_reciprocal_rss<float,float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_reciprocal_rss<float,float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_reciprocal_rss<float_complext, float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_correlation<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_correlation<float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_axpy<float>( cuNDArray<float>*, cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<float,1>( uintd1::Type, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<complext<float>,1>( uintd1::Type, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<complext<float>,2>( uintd2::Type, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<complext<float>,3>( uintd3::Type, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<complext<float>,4>( uintd4::Type, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,1>,1>( uintd1::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,2>,1>( uintd1::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,3>,1>( uintd1::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,4>,1>( uintd1::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,1>,2>( uintd2::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,2>,2>( uintd2::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,3>,2>( uintd2::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,4>,2>( uintd2::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,1>,3>( uintd3::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,2>,3>( uintd3::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,3>,3>( uintd3::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,4>,3>( uintd3::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,1>,4>( uintd4::Type, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,2>,4>( uintd4::Type, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,3>,4>( uintd4::Type, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<float,4>,4>( uintd4::Type, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float,1>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float_complext,1>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float,2>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float_complext,2>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float,3>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float_complext,3>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float,4>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<float_complext,4>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
cuNDA_real_to_complext<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_complext_to_real<float>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_downsample<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_downsample<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_downsample<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_downsample<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_nn<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_nn<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_nn<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_nn<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_lin<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_lin<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_lin<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > 
cuNDA_upsample_lin<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_clear<float>( cuNDArray<float>*,float, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_clear<float_complext>( cuNDArray<float_complext>*,float_complext, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_reciprocal<float>( cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_reciprocal<float_complext>( cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_sqrt<float>( cuNDArray<float>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_reciprocal_sqrt<float>( cuNDArray<float>*, cuNDA_device );



template EXPORTGPUCORE bool cuNDA_abs<float>( cuNDArray<float>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_abs<floatd1::Type>( cuNDArray<floatd1::Type>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_abs<floatd2::Type>( cuNDArray<floatd2::Type>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_abs<floatd3::Type>( cuNDArray<floatd3::Type>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_abs<floatd4::Type>( cuNDArray<floatd4::Type>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_threshold_min<float>(float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_threshold_max<float>(float, cuNDArray<float>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_rss_normalize<float,float>( cuNDArray<float>*, unsigned int, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_rss_normalize<float,float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device );


template EXPORTGPUCORE bool cuNDA_add<float_complext>( float_complext, cuNDArray<float_complext>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_add<float>( float, cuNDArray<float>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scale<float>( float, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scale<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_scale<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_scale_conj<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scale<float>( cuNDArray<float>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_axpy<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_axpy<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,1>(uintd1::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float_complext,1>(uintd1::Type, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,2>(uintd2::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float_complext,2>(uintd2::Type, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,3>(uintd3::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float_complext,3>(uintd3::Type, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,4>(uintd4::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float_complext,4>(uintd4::Type, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float,1>(floatd1::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float_complext,1>(floatd1::Type, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float,2>(floatd2::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float_complext,2>(floatd2::Type, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float,3>(floatd3::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float_complext,3>(floatd3::Type, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float,4>(floatd4::Type, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<float,float_complext,4>(floatd4::Type, cuNDArray<float_complext>*, cuNDA_device );


template EXPORTGPUCORE float cuNDA_dot<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE float_complext cuNDA_dot<float_complext>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE float cuNDA_asum<float,float>( cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE float cuNDA_asum<float,float_complext>( cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_axpy<float>( float, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_axpy<float_complext>( float_complext, cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scal<float>( float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_scal<float_complext>( float_complext, cuNDArray<float_complext>*, cuNDA_device );



template EXPORTGPUCORE bool cuNDA_shrink1<float,float>( float, cuNDArray<float>*, cuNDArray<float>* );
template EXPORTGPUCORE bool cuNDA_shrink1<float,float_complext>( float, cuNDArray<float_complext>*, cuNDArray<float_complext>* );

template EXPORTGPUCORE bool cuNDA_shrinkd<float,float>( float, cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>* );
template EXPORTGPUCORE bool cuNDA_shrinkd<float,float_complext>( float, cuNDArray<float>*, cuNDArray<float_complext>*, cuNDArray<float_complext>* );

template EXPORTGPUCORE 
bool cuNDA_origin_mirror<float,1>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
bool cuNDA_origin_mirror<float,2>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
bool cuNDA_origin_mirror<float,3>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE 
bool cuNDA_origin_mirror<float,4>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);

template EXPORTGPUCORE bool
cuNDA_origin_mirror<float_complext,1>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE bool 
cuNDA_origin_mirror<float_complext,2>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE bool 
cuNDA_origin_mirror<float_complext,3>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE bool 
cuNDA_origin_mirror<float_complext,4>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);


// Instanciation -- double precision

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_sum<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<complext< double> > >
cuNDA_sum<complext<double> >( cuNDArray<complext< double> >*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<1>::Type> > 
cuNDA_sum<doubled<1>::Type>( cuNDArray<doubled<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<2>::Type> > 
cuNDA_sum<doubled<2>::Type>( cuNDArray<doubled<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<3>::Type> >
 cuNDA_sum<doubled<3>::Type>( cuNDArray<doubled<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<4>::Type> >
 cuNDA_sum<doubled<4>::Type>( cuNDArray<doubled<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_expand<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_expand<double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device);


template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_cAbs<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
cuNDA_cAbs<double,double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
cuNDA_cNorm<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
cuNDA_cNorm<double,double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm<double,1>( cuNDArray<doubled<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm<double,2>( cuNDArray<doubled<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm<double,3>( cuNDArray<doubled<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm<double,4>( cuNDArray<doubled<4>::Type>*, cuNDA_device, cuNDA_device );
/*
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm_squared<double,double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
 cuNDA_norm_squared<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );
*/
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm_squared<double,1>( cuNDArray<doubled<1>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm_squared<double,2>( cuNDArray<doubled<2>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm_squared<double,3>( cuNDArray<doubled<3>::Type>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_norm_squared<double,4>( cuNDArray<doubled<4>::Type>*, cuNDA_device, cuNDA_device );


template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_ss<double,double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_ss<double,double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_ss<double_complext, double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_rss<double,double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_rss<double,double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_rss<double_complext, double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_reciprocal_rss<double,double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_reciprocal_rss<double,double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_reciprocal_rss<double_complext, double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_correlation<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template<> EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_correlation<double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_axpy<double>( cuNDArray<double>*, cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<double,1>( uintd1::Type, cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<complext<double> ,1>( uintd1::Type, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<complext<double> ,2>( uintd2::Type, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<complext<double> ,3>( uintd3::Type, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<complext<double> ,4>( uintd4::Type, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );



template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,1>,1>( uintd1::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,2>,1>( uintd1::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,3>,1>( uintd1::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,4>,1>( uintd1::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,1>,2>( uintd2::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,2>,2>( uintd2::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,3>,2>( uintd2::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,4>,2>( uintd2::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,1>,3>( uintd3::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,2>,3>( uintd3::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,3>,3>( uintd3::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,4>,3>( uintd3::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,1>,4>( uintd4::Type, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,2>,4>( uintd4::Type, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,3>,4>( uintd4::Type, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_crop<vector_td<double,4>,4>( uintd4::Type, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double,1>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double_complext,1>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double,2>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double_complext,2>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double,3>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double_complext,3>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double,4>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool
cuNDA_expand_with_zero_fill<double_complext,4>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
cuNDA_real_to_complext<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_complext_to_real<double>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_downsample<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_downsample<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_downsample<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_downsample<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_nn<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_nn<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_nn<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_nn<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_lin<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_lin<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_lin<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > 
cuNDA_upsample_lin<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_clear<double>( cuNDArray<double>*,double, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_clear<double_complext>( cuNDArray<double_complext>*,double_complext, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_reciprocal<double>( cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_reciprocal<double_complext>( cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_sqrt<double>( cuNDArray<double>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_reciprocal_sqrt<double>( cuNDArray<double>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_abs<double>( cuNDArray<double>*, cuNDA_device );


template EXPORTGPUCORE bool cuNDA_abs<doubled1::Type>( cuNDArray<doubled1::Type>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_abs<doubled2::Type>( cuNDArray<doubled2::Type>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_abs<doubled3::Type>( cuNDArray<doubled3::Type>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_abs<doubled4::Type>( cuNDArray<doubled4::Type>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_threshold_min<double>(double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_threshold_max<double>(double, cuNDArray<double>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_rss_normalize<double,double>( cuNDArray<double>*, unsigned int, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_rss_normalize<double,double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device );


template EXPORTGPUCORE bool cuNDA_add<double_complext>( double_complext, cuNDArray<double_complext>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_add<double>( double, cuNDArray<double>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scale<double>( double, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scale<double>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_scale<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scale<double>( cuNDArray<double>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scale_conj<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_axpy<double>( cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_axpy<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,1>(uintd1::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double_complext,1>(uintd1::Type, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,2>(uintd2::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double_complext,2>(uintd2::Type, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,3>(uintd3::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double_complext,3>(uintd3::Type, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,4>(uintd4::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double_complext,4>(uintd4::Type, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double,1>(doubled1::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double_complext,1>(doubled1::Type, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double,2>(doubled2::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double_complext,2>(doubled2::Type, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double,3>(doubled3::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double_complext,3>(doubled3::Type, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double,4>(doubled4::Type, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_zero_fill_border<double,double_complext,4>(doubled4::Type, cuNDArray<double_complext>*, cuNDA_device );


template EXPORTGPUCORE double cuNDA_dot<double>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE double_complext cuNDA_dot<double_complext>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE double cuNDA_asum<double,double>( cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE double cuNDA_asum<double,double_complext>( cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_axpy<double>( double, cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_axpy<double_complext>( double_complext, cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_scal<double>( double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE bool cuNDA_scal<double_complext>( double_complext, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE bool cuNDA_shrink1<double,double>( double, cuNDArray<double>*, cuNDArray<double>* );
template EXPORTGPUCORE bool cuNDA_shrink1<double,double_complext>( double, cuNDArray<double_complext>*, cuNDArray<double_complext>* );

template EXPORTGPUCORE bool cuNDA_shrinkd<double,double>( double, cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>* );
template EXPORTGPUCORE bool cuNDA_shrinkd<double,double_complext>( double, cuNDArray<double>*, cuNDArray<double_complext>*, cuNDArray<double_complext>* );

template EXPORTGPUCORE 
bool cuNDA_origin_mirror<double,1>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
bool cuNDA_origin_mirror<double,2>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
bool cuNDA_origin_mirror<double,3>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE 
bool cuNDA_origin_mirror<double,4>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);

template EXPORTGPUCORE bool
cuNDA_origin_mirror<double_complext,1>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE bool 
cuNDA_origin_mirror<double_complext,2>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE bool 
cuNDA_origin_mirror<double_complext,3>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE bool 
cuNDA_origin_mirror<double_complext,4>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
