/*
  CUDA implementation of the NFFT.

  -----------

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12): 1974-1985. 
*/

// Includes - our own code
#include "NFFT.h"
#include "cuNDFFT.h"

#include "cuNDArray.h"
#include "ndarray_vector_td_utilities.h"
#include "vector_td_utilities.h"

#include "check_CUDA.h"

// Includes - CUDA
#include <math_constants.h>
#include <cufft.h>

// Includes - CUDA thrust
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

// Includes - stdlibs
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

using namespace std;
using namespace thrust;

//
// Some defines used to configure our kernels
//

// TODO: make "d x generation"-dimensional struct to contain MAX_COILS and THREADS_PER_KERNEL for the NFFT and NFFT_H respectively

// Block configuration for kernel
// ------------------------------

#define NFFT_MAX_COILS_2D_COMPUTE_1x	                                          8
#define NFFT_MAX_COILS_2D_COMPUTE_2x	                                         16
#define NFFT_THREADS_PER_2D_KERNEL_1x						256
#define NFFT_THREADS_PER_2D_KERNEL_2x						256

#define NFFT_H_MAX_COILS_2D_COMPUTE_1x	                                          8
#define NFFT_H_MAX_COILS_2D_COMPUTE_2x	                                         16
#define NFFT_H_THREADS_PER_2D_KERNEL_1x						256
#define NFFT_H_THREADS_PER_2D_KERNEL_2x						256

#define NFFT_MAX_COILS_3D_COMPUTE_1x	                                          8
#define NFFT_MAX_COILS_3D_COMPUTE_2x	                                         16
#define NFFT_THREADS_PER_3D_KERNEL_1x						256
#define NFFT_THREADS_PER_3D_KERNEL_2x						256

#define NFFT_H_MAX_COILS_3D_COMPUTE_1x	                                          8
#define NFFT_H_MAX_COILS_3D_COMPUTE_2x	                                         16
#define NFFT_H_THREADS_PER_3D_KERNEL_1x						256
#define NFFT_H_THREADS_PER_3D_KERNEL_2x						256

#define NFFT_MAX_COILS_4D_COMPUTE_1x	                                          8
#define NFFT_MAX_COILS_4D_COMPUTE_2x	                                         16
#define NFFT_THREADS_PER_4D_KERNEL_1x						256
#define NFFT_THREADS_PER_4D_KERNEL_2x						256

#define NFFT_H_MAX_COILS_4D_COMPUTE_1x	                                          8
#define NFFT_H_MAX_COILS_4D_COMPUTE_2x	                                         16
#define NFFT_H_THREADS_PER_4D_KERNEL_1x						256
#define NFFT_H_THREADS_PER_4D_KERNEL_2x						256
  
//
// Reference to shared memory
//

extern __shared__ char _shared_mem[];

// 
// Includes containing the NFFT convolution implementation
// 

#include "KaiserBessel_kernel.cu"
#include "NFFT_kernel.cu"
#include "NFFT_H_kernel.cu"
#include "NFFT_preprocess_kernel.cu"

// Some device properties we query once to eliminate runtime overhead
//
static int num_devices = 0;
static int *warp_size = 0x0;
static int *major = 0x0;

// Default template arguments seems to require c++-0x. 
typedef float dummy;

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
  major = new int[num_devices];

  if( !warp_size || !major ) {
    cout << endl << "Error: trivial malloc failed!" << endl ;
    return false;
  }

  for( int device=0; device<num_devices; device++ ){

    cudaDeviceProp deviceProp; 
    
    if( cudaGetDeviceProperties( &deviceProp, device ) != cudaSuccess) {
      cout << endl << "Error: unable to determine device properties." << endl ;
      return false;
    }

    warp_size[device] = deviceProp.warpSize;
    major[device] = deviceProp.major;
  }

  if( cudaSetDevice(old_device) != cudaSuccess ) {
    cerr << endl << "Error: unable to restore device no";
    return false;
  }

  return true;
}


// Common multi-device handling: prepare
//
template<class I1, class I2, class I3>
static bool prepare( int device, int *old_device, 
		     cuNDArray<I1> *in1,       cuNDArray<I1> **in1_int,
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> **in2_int = 0x0,
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> **in3_int = 0x0 )
{
  // Get current Cuda device
  if( cudaGetDevice(old_device) != cudaSuccess ) {
    cerr << endl << "unable to get device no";
    return false;
  }

  if( device != *old_device && cudaSetDevice(device) != cudaSuccess) {
    cerr << endl << "unable to set device no";
    return false;
  }
  
  // Transfer arrays to compute device if necessary
  if( in1 ){
    if( device != in1->get_device() )
      *in1_int = new cuNDArray<I1>(*in1); // device transfer
    else
      *in1_int = in1;
  }
  
  if( in2 ){
    if( device != in2->get_device() )
      *in2_int = new cuNDArray<I2>(*in2); // device transfer
    else
      *in2_int = in2;
  }

  if( in3 ){
    if( device != in3->get_device() )
      *in3_int = new cuNDArray<I3>(*in3); // device transfer
    else
      *in3_int = in3;
  }
  
  return true;
}  

// Common multi-device handling: restore
//
template<class I1, class I2, class I3>
static bool restore( int old_device, cuNDArray<I1> *out, 
		     cuNDArray<I1> *in1, cuNDArray<I1> *in1_int,
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> *in2_int = 0x0,
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> *in3_int = 0x0 )
{
  if( in1 && out && out->get_device() != in1_int->get_device() ){ 
    *out = *in1_int; } // device transfer by assignment
  
  // Check if internal array needs deletion (they do only if they were created in ::prepare()
  //
  if( in1 && in1->get_device() != in1_int->get_device() ){
    delete in1_int;
  }   
  if( in2 && in2->get_device() != in2_int->get_device() ){
    delete in2_int;
  }   
  if( in3 && in3->get_device() != in3_int->get_device() ){
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

//
// Public class methods
//

template<class REAL, unsigned int D> 
NFFT_plan<REAL,D>::NFFT_plan()
{
  // Minimal initialization
  barebones();
}

template<class REAL, unsigned int D> 
NFFT_plan<REAL,D>::NFFT_plan( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W, int device )
{
  // Minimal initialization
  barebones();

  // Setup plan
  setup( matrix_size, matrix_size_os, W, device );
}

template<class REAL, unsigned int D> 
NFFT_plan<REAL,D>::~NFFT_plan()
{
  wipe(NFFT_WIPE_ALL);
}

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W, int _device )
{
  // Free memory
  wipe(NFFT_WIPE_ALL);

  //
  // Check if the device is valid
  //

  if( _device<0 ){
    if( cudaGetDevice( &device ) != cudaSuccess ){
      cout << endl << "NFFT_plan::setup: unable to determine device properties." << endl ;
      return false;
    }
  }
  else
    device = _device;
  
  if( !initialize_static_variables() ){
    cout << endl << "NFFT_plan::setup: unable to query device properties." << endl ;
    return false;
  }

  typename uintd<D>::Type vec_warp_size = to_vector_td<unsigned int,D>(warp_size[device]);

  //
  // Check input against certain requirements
  //
  
  if( sum(matrix_size%vec_warp_size) || sum(matrix_size_os%vec_warp_size) ){
    cout << endl << "Illegal matrix size for the NFFT plan (not a multiple of the warp size)" << endl;
    return false;
  }

  //
  // Setup private variables
  //

  this->matrix_size = matrix_size;
  this->matrix_size_os = matrix_size_os;

  REAL W_half = get_half<REAL>()*W;
  vector_td<REAL,D> W_vec = to_vector_td<REAL,D>(W_half);

  matrix_size_wrap = to_uintd<REAL,D>( ceil(W_vec));
  matrix_size_wrap<<=1; 
  
  alpha = (REAL) matrix_size_os.vec[0] / (REAL) matrix_size.vec[0];
  
  REAL one = get_one<REAL>();
  if( alpha < one ){
    cout << endl << "Illegal oversampling ratio suggested" << endl;
    return false;
  }

  this->W = W;

  unsigned int fracX = matrix_size_os.vec[0] / matrix_size.vec[0];
  unsigned int moduX = matrix_size_os.vec[0] % matrix_size.vec[0];

  for( unsigned int dim=0; dim<D; dim++){
    
    if( ((matrix_size_os.vec[dim]/(matrix_size.vec[dim])) != fracX || (matrix_size_os.vec[dim]%(matrix_size.vec[dim]) != moduX) )){
      cout << endl << "Oversampling ratio is not constant between dimensions" << endl;
      return false;
    }
  }
  
  // Compute Kaiser-Bessel beta
  bool success = 
    compute_beta();
  
  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to get device no" << endl;
  }  
  if( device != device_no_old && cudaSetDevice(device) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }  

  initialized = success;  

  // Calculate deapodization filter
  if( success ) success = compute_deapodization_filter();
  
  initialized = success;  

  if( device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to restore device" << endl;
  }
  
  return success;
}

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory, NFFT_prep_mode mode )
{
  if( !trajectory || trajectory->get_number_of_elements()==0 ){
    cout << endl << "NFFT_plan::preprocess: 0x0 or empty trajectory provided." << endl;
    return false;
  }
  
  if( !initialized ){
    cout << endl << "NFFT_plan::preprocess: NFFT_plan::setup must be invoked prior to preprocessing." << endl;
    return false;
  }

  if( !wipe(NFFT_WIPE_PREPROCESSING) ){
    cout << endl << "NFFT_plan::preprocess: unable to clean plan." << endl;
    return false;
  }

  cuNDArray<typename reald<REAL,D>::Type> *trajectory_int;
  int old_device;

  if( !prepare<typename reald<REAL,D>::Type,dummy,dummy>(device, &old_device, trajectory, &trajectory_int ) ){
    cout << endl << "NFFT_plan::preprocess: device preparation error." << endl;
    return false;
  }
    
  number_of_samples = trajectory_int->get_size(0);
  number_of_frames = trajectory_int->get_number_of_elements()/number_of_samples;
  
  // Make Thrust device vector of trajectory and samples
  device_vector< vector_td<REAL,D> > trajectory_positions_in( device_pointer_cast< vector_td<REAL,D> >(trajectory_int->get_data_ptr()), 
							      device_pointer_cast< vector_td<REAL,D> >(trajectory_int->get_data_ptr()+trajectory_int->get_number_of_elements() ));

  trajectory_positions = new device_vector< vector_td<REAL,D> >( trajectory_int->get_number_of_elements() );
    
  CHECK_FOR_CUDA_ERROR();

  vector_td<REAL,D> matrix_size_os_real = to_reald<REAL,unsigned int,D>( matrix_size_os );
  vector_td<REAL,D> matrix_size_os_plus_wrap_real = to_reald<REAL,unsigned int,D>( (matrix_size_os+matrix_size_wrap)>>1 );

  // convert input trajectory in [-1/2;1/2] to [0;matrix_size_os]
  thrust::transform( trajectory_positions_in.begin(), trajectory_positions_in.end(), trajectory_positions->begin(), trajectory_scale<REAL,D>(matrix_size_os_real, matrix_size_os_plus_wrap_real) );
  
  if( mode != NFFT_PREP_FORWARDS ){
    
    // allocate storage for and compute temporary prefix-sum variable (#cells influenced per sample)
    device_vector<unsigned int> c_p_s(trajectory_int->get_number_of_elements());
    device_vector<unsigned int> c_p_s_ps(trajectory_int->get_number_of_elements());
    
    CHECK_FOR_CUDA_ERROR();
    
    REAL half_W = get_half<REAL>()*W;
    thrust::plus<unsigned int> binary_op;
    thrust::transform(trajectory_positions->begin(), trajectory_positions->end(), c_p_s.begin(), compute_num_cells_per_sample<REAL,D>(half_W));
    inclusive_scan( c_p_s.begin(), c_p_s.end(), c_p_s_ps.begin(), binary_op ); // prefix sum
    
    // Build the vector of (grid_idx, sample_idx) tuples. Actually kept in two seperate vectors.
    unsigned int num_pairs = c_p_s_ps.back();
    c_p_s.clear();

    thrust::device_vector<unsigned int> *tuples_first = new device_vector<unsigned int>(num_pairs);
    tuples_last = new device_vector<unsigned int>(num_pairs);
    
    CHECK_FOR_CUDA_ERROR();
    
    // Fill tuple vector
    write_pairs<REAL,D>( matrix_size_os, matrix_size_wrap, number_of_samples, number_of_frames, W, 
			 raw_pointer_cast(&(*trajectory_positions)[0]), raw_pointer_cast(&c_p_s_ps[0]), 
			 raw_pointer_cast(&(*tuples_first)[0]), raw_pointer_cast(&(*tuples_last)[0]) );
    c_p_s_ps.clear();

    // Sort by grid indices
    sort_by_key( tuples_first->begin(), tuples_first->end(), tuples_last->begin() );
    
    // each bucket_begin[i] indexes the first element of bucket i's list of points
    // each bucket_end[i] indexes one past the last element of bucket i's list of points

    bucket_begin = new device_vector<unsigned int>(number_of_frames*prod(matrix_size_os+matrix_size_wrap));
    bucket_end   = new device_vector<unsigned int>(number_of_frames*prod(matrix_size_os+matrix_size_wrap));
    
    CHECK_FOR_CUDA_ERROR();
    
    // find the beginning of each bucket's list of points
    counting_iterator<unsigned int> search_begin(0);
    lower_bound(tuples_first->begin(), tuples_first->end(), search_begin, search_begin + number_of_frames*prod(matrix_size_os+matrix_size_wrap), bucket_begin->begin() );
    
    // find the end of each bucket's list of points
    upper_bound(tuples_first->begin(), tuples_first->end(), search_begin, search_begin + number_of_frames*prod(matrix_size_os+matrix_size_wrap), bucket_end->begin() );
  
    delete tuples_first;
  }
    
  preprocessed_NFFT = true;

  if( mode != NFFT_PREP_FORWARDS )
    preprocessed_NFFT_H = true;

  if( !restore<typename reald<REAL,D>::Type,dummy,dummy>(old_device, trajectory, trajectory, trajectory_int) ){
    cout << endl << "NFFT_plan::preprocess: unable to restore compute device." << endl;
    return false;
  }
  
  return true;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::compute( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, 
			    cuNDArray<REAL> *weights, NFFT_comp_mode mode )
{  

  // Validity checks
  
  unsigned char components;

  if( mode == NFFT_FORWARDS ) 
    components = NFFT_CONVOLUTION + NFFT_FFT + NFFT_DEAPODIZATION;

  else // backwards
    components = NFFT_H_CONVOLUTION + NFFT_FFT + NFFT_DEAPODIZATION;

  if( !check_consistency( samples, image, weights, components ) ){
    cout << endl << "NFFT_plan::compute: input consistency check failed." << endl;
    return false;
  }
  
  cuNDArray<typename complext<REAL>::Type> *samples_int = 0x0, *image_int = 0x0;
  cuNDArray<REAL> *weights_int = 0x0;
  int old_device;

  if( !prepare<typename complext<REAL>::Type,typename complext<REAL>::Type,REAL>
      (device, &old_device, samples, &samples_int, image, &image_int, weights, &weights_int ) ){
    cout << endl << "NFFT_plan::compute: device preparation error." << endl;
    return false;
  }

  bool success;  
  cuNDArray<typename complext<REAL>::Type> *working_image = 0x0, *working_samples = 0x0;
  
  typename uintd<D>::Type image_dims = vector_to_uintd<D>(*image->get_dimensions()); 
  bool oversampled_image = (image_dims==matrix_size_os); 

  vector<unsigned int> vec_dims = uintd_to_vector<D>(matrix_size_os); 
  for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
    vec_dims.push_back(image->get_size(d));
  
  switch(mode){

  case NFFT_FORWARDS:
    
    if( !oversampled_image ){
      working_image = new cuNDArray<typename complext<REAL>::Type>; working_image->create(&vec_dims);
      cuNDA_expand_with_zero_fill<typename complext<REAL>::Type,D>( image_int, working_image );
    }
    else{
      working_image = image_int;
    }
    
    success = compute_NFFT( samples_int, working_image );

    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    break;
    
  case NFFT_BACKWARDS:

    // Density compensation
    if( weights_int ){
      working_samples = new cuNDArray<typename complext<REAL>::Type>(*samples_int);    
      cuNDA_scale( weights_int, working_samples );
    }
    else{
      working_samples = samples_int;
    }
    
    if( !oversampled_image ){
      working_image = new cuNDArray<typename complext<REAL>::Type>; working_image->create(&vec_dims); 
    }
    else{
      working_image = image_int;
    }

    success = compute_NFFT_H( working_samples, working_image );

    if( success ){
      if( !oversampled_image ){
	cuNDA_crop<typename complext<REAL>::Type,D>( (matrix_size_os-matrix_size)>>1, working_image, image_int );
      }
    }
    
    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    if( weights ){
      delete working_samples; working_samples = 0x0;
    }
    
    break;
  };

  if( !restore<typename complext<REAL>::Type,typename complext<REAL>::Type,REAL>
      (old_device, (mode==NFFT_FORWARDS)?samples:image, 
       (mode==NFFT_FORWARDS)?samples:image, (mode==NFFT_FORWARDS)?samples_int:image_int,
       (mode!=NFFT_FORWARDS)?samples:image, (mode!=NFFT_FORWARDS)?samples_int:image_int,
       weights, weights_int ) ){
    cout << endl << "NFFT_plan::compute: unable to restore compute device." << endl;
    return false;
  }
  
  CHECK_FOR_CUDA_ERROR();
  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::compute_iteration( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, 
				      cuNDArray<REAL> *weights, NFFT_comp_mode mode )
{  
  // Validity checks
  
  unsigned char components = NFFT_CONVOLUTION + NFFT_H_CONVOLUTION + NFFT_FFT + NFFT_DEAPODIZATION;
  
  if( !check_consistency( samples, image, weights, components ) )
    return false;
  
  cuNDArray<typename complext<REAL>::Type> *samples_int = 0x0, *image_int = 0x0;
  cuNDArray<REAL> *weights_int = 0x0;
  int old_device;

  if( !prepare<typename complext<REAL>::Type,typename complext<REAL>::Type,REAL>
      (device, &old_device, samples, &samples_int, image, &image_int, weights, &weights_int ) ){
    cout << endl << "NFFT_plan::compute_iteration: device preparation error." << endl;
    return false;
  }

  bool success;  
  cuNDArray<typename complext<REAL>::Type> *working_image = 0x0, *working_samples = 0x0;

  typename uintd<D>::Type image_dims = vector_to_uintd<D>(*image->get_dimensions()); 
  bool oversampled_image = (image_dims==matrix_size_os); 
 
  vector<unsigned int> vec_dims = uintd_to_vector<D>(matrix_size_os); 
  for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
    vec_dims.push_back(image->get_size(d));
  
  switch(mode){

  case NFFT_FORWARDS: // iteration from image
    
    if( !oversampled_image ){
      working_image = new cuNDArray<typename complext<REAL>::Type>; working_image->create(&vec_dims); 
      cuNDA_expand_with_zero_fill<typename complext<REAL>::Type,D>( image_int, working_image );
    }
    else{
      working_image = image_int;
    }
    
    working_samples = samples;

    success = compute_NFFT( working_samples, working_image );
    
    if( success ){
      
      // Density compensation
      if( weights ){
	cuNDA_scale( weights_int, working_samples );
      }
      
      success = compute_NFFT_H( working_samples, working_image ); 
    }      
    
    if( success ){
      if( !oversampled_image ){
	cuNDA_crop<typename complext<REAL>::Type,D>( (matrix_size_os-matrix_size)>>1, working_image, image_int );
	delete working_image; working_image = 0x0;
      }
    }
    
    break;

  case NFFT_BACKWARDS: // iteration from samples

    // Density compensation
    if( weights ){
      working_samples = new cuNDArray<typename complext<REAL>::Type>(*samples_int);
      cuNDA_scale( weights_int, working_samples );
    }
    else{
      working_samples = samples_int;
    }
    
    if( !oversampled_image ){  
      working_image = new cuNDArray<typename complext<REAL>::Type>; working_image->create(&vec_dims); 
    } 
    else
      working_image = image_int;

    success = compute_NFFT_H( working_samples, working_image );
    
    if( success ){
      cuNDA_zero_fill_border<typename complext<REAL>::Type,D>( matrix_size, working_image );      
      success = compute_NFFT( samples_int, working_image );
    }
  
    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    if( weights ){
      delete working_samples; working_samples = 0x0;
    }
    
    break;
  };

  if( !restore<typename complext<REAL>::Type,typename complext<REAL>::Type,REAL>
      (old_device, (mode==NFFT_FORWARDS)?samples:image, 
       (mode==NFFT_FORWARDS)?samples:image, (mode==NFFT_FORWARDS)?samples_int:image_int,
       (mode!=NFFT_FORWARDS)?samples:image, (mode!=NFFT_FORWARDS)?samples_int:image_int,
       weights, weights_int ) ){
    cout << endl << "NFFT_plan::compute_iteration: unable to restore compute device." << endl;
    return false;
  }

  CHECK_FOR_CUDA_ERROR();
  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::convolve( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, 
			     cuNDArray<REAL> *weights, NFFT_comp_mode mode, bool accumulate )
{
  unsigned char components;

  if( mode == NFFT_FORWARDS ) 
    components = NFFT_CONVOLUTION;
  else
    components = NFFT_H_CONVOLUTION;
  
  if( !check_consistency( samples, image, weights, components ) )
    return false;
  
  cuNDArray<typename complext<REAL>::Type> *samples_int = 0x0, *image_int = 0x0;
  cuNDArray<REAL> *weights_int = 0x0;
  int old_device;

  if( !prepare<typename complext<REAL>::Type,typename complext<REAL>::Type,REAL>
      (device, &old_device, samples, &samples_int, image, &image_int, weights, &weights_int ) ){
    cout << endl << "NFFT_plan::convolve: device preparation error." << endl;
    return false;
  }

  bool success;  
  cuNDArray<typename complext<REAL>::Type> *working_samples = 0x0;
  
  typename uintd<D>::Type image_dims = vector_to_uintd<D>(*image->get_dimensions()); 
  bool oversampled_image = (image_dims==matrix_size_os); 
 
  if( !oversampled_image ){
    cout << endl << "NFFT_plan::convolve: ERROR: oversampled image not provided as input." << endl;
    return false;
  }

  vector<unsigned int> vec_dims = uintd_to_vector<D>(matrix_size_os); 
  for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
    vec_dims.push_back(image->get_size(d));
  
  switch(mode){

  case NFFT_FORWARDS:
    success = convolve_NFFT( samples_int, image_int, accumulate );
    break;
    
  case NFFT_BACKWARDS:

    // Density compensation
    if( weights ){
      working_samples = new cuNDArray<typename complext<REAL>::Type>(*samples_int);
      cuNDA_scale( weights_int, working_samples );
    }
    else{
      working_samples = samples_int;
    }
    
    success = convolve_NFFT_H( working_samples, image_int, accumulate );

    if( weights ){
      delete working_samples; working_samples = 0x0;
    }
    
    break;
  }

  if( !restore<typename complext<REAL>::Type,typename complext<REAL>::Type,REAL>
      (old_device, (mode==NFFT_FORWARDS)?samples:image, 
       (mode==NFFT_FORWARDS)?samples:image, (mode==NFFT_FORWARDS)?samples_int:image_int,
       (mode!=NFFT_FORWARDS)?samples:image, (mode!=NFFT_FORWARDS)?samples_int:image_int,
       weights, weights_int ) ){
    cout << endl << "NFFT_plan::convolve: unable to restore compute device." << endl;
    return false;
  }

  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::fft(cuNDArray<typename complext<REAL>::Type> *data, NFFT_comp_mode mode, bool do_scale )
{
  cuNDArray<typename complext<REAL>::Type> *data_int = 0x0;
  int old_device;
  
  if( !prepare<typename complext<REAL>::Type,dummy,dummy>(device, &old_device, data, &data_int ) ){
    cout << endl << "NFFT_plan::fft: device preparation error." << endl;
    return false;
  }
  
  typename uintd<D>::Type _dims_to_transform = counting_vec<D>();
  vector<unsigned int> dims_to_transform = uintd_to_vector<D>( _dims_to_transform );
  
  int res;
  if( mode == NFFT_FORWARDS ){
    res = cuNDFFT<typename complext<REAL>::Type>().fft( data_int, &dims_to_transform ); 
  }
  else{
    res = cuNDFFT<typename complext<REAL>::Type>().ifft( data_int, &dims_to_transform, do_scale );
  }

  if( !restore<typename complext<REAL>::Type,dummy,dummy>(old_device, data, data, data_int) ){
    cout << endl << "NFFT_plan::fft: unable to restore compute device." << endl;
    return false;
  }

  if( res == 0 )
    return true;
  else
    return false;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::deapodize( cuNDArray<typename complext<REAL>::Type> *image )
{
  unsigned char components;
  
  components = NFFT_FFT;
  
  if( !check_consistency( 0x0, image, 0x0, components ) )
    return false;

  cuNDArray<typename complext<REAL>::Type> *image_int = 0x0;
  int old_device;
  
  if( !prepare<typename complext<REAL>::Type,dummy,dummy>(device, &old_device, image, &image_int ) ){
    cout << endl << "NFFT_plan::deapodize: device preparation error." << endl;
    return false;
  }

  typename uintd<D>::Type image_dims = vector_to_uintd<D>(*image->get_dimensions()); 
  bool oversampled_image = (image_dims==matrix_size_os); 
  
  if( !oversampled_image ){
    cout << endl << "NFFT_plan::deapodize: ERROR: oversampled image not provided as input." << endl;
    return false;
  }
  
  cuNDA_scale( deapodization_filter.get(), image_int );
  
  if( !restore<typename complext<REAL>::Type,dummy,dummy>(old_device, image, image, image_int) ){
    cout << endl << "NFFT_plan::deapodize: unable to restore compute device." << endl;
    return false;
  }

  return true;
}

//
// Private class methods
//

template<class REAL, unsigned int D> bool 
NFFT_plan<REAL,D>::check_consistency( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, unsigned char components )
{
  if( !initialized ){
    cout << endl << "NFFT_plan: Unable to proceed without setup." << endl;
    return false;
  }
  
  if( (components & NFFT_CONVOLUTION ) && !preprocessed_NFFT ){
    cout << endl << "NFFT_plan: Unable to compute NFFT before preprocessing." << endl;
    return false;
  }
  
  if( (components & NFFT_H_CONVOLUTION ) && !preprocessed_NFFT_H ){
    cout << endl << "NFFT_plan: Unable to compute NFFT before preprocessing." << endl;
    return false;
  }
  
  if( ((components & NFFT_CONVOLUTION ) || (components & NFFT_H_CONVOLUTION )) && !(image && samples) ){
    cout << endl << "NFFT_plan: Unable to process 0x0 input/output." << endl;
    return false;
  }
  
  if( ((components & NFFT_FFT) || (components & NFFT_DEAPODIZATION )) && !image ){
    cout << endl << "NFFT_plan: Unable to process 0x0 input." << endl;
    return false;
  }

  if( image->get_number_of_dimensions() < D ){
    cout << endl << "NFFT_plan: Number of image dimensions mismatch the plan." << endl;
    return false;
  }    

  typename uintd<D>::Type image_dims = vector_to_uintd<D>( *image->get_dimensions() );
  bool oversampled_image = (image_dims==matrix_size_os);
  
  if( !((oversampled_image) ? (image_dims == matrix_size_os) : (image_dims == matrix_size) )){
    cout << endl << "NFFT_plan: Image dimensions mismatch." << endl;
    return false;
  }
  
  if( (components & NFFT_CONVOLUTION ) || (components & NFFT_H_CONVOLUTION ) ){
    
    if( (samples->get_number_of_elements() == 0) || (samples->get_number_of_elements() % (number_of_frames*number_of_samples)) ){
      cout << endl << "NFFT_plan: The number of samples is not a multiple of #samples/frame x #frames as requested through preprocessing" << endl;
      return false;
    }
    
    unsigned int num_batches_in_samples_array = samples->get_number_of_elements()/(number_of_frames*number_of_samples);

    unsigned int num_batches_in_image_array = 1;
    for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ ){
      num_batches_in_image_array *= image->get_size(d);
    }
    num_batches_in_image_array /= number_of_frames;

    if( num_batches_in_samples_array != num_batches_in_image_array ){
      cout << endl << "NFFT_plan: Number of batches mismatch between samples and image arrays" << endl;
      return false;
    }
  }
  
  if( components & NFFT_H_CONVOLUTION ){
    if( weights ){
 
      if( weights->get_number_of_elements() == 0 ||
	  !( weights->get_number_of_elements() == number_of_samples || 
	     weights->get_number_of_elements() == number_of_frames*number_of_samples) ){

	cout << endl << "NFFT_plan: The number of weights should match #samples/frame x #frames as requested through preprocessing" << endl;
	return false;
      }
    }
  }
  
  return true;
}

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::barebones()
{	
  // These are the fundamental booleans checked before accessing the various member pointers
  initialized = preprocessed_NFFT = preprocessed_NFFT_H = false;
  
  // and specify the device
  if (cudaGetDevice(&device) != cudaSuccess) {
    cerr << "NFFT_plan::barebones:: unable to get device no" << endl;
    return false;
  }

  return true;
}

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::wipe( NFFT_wipe_mode mode )
{
  // Get current Cuda device
  int old_device;
  if( cudaGetDevice(&old_device) != cudaSuccess ) {
    cerr << endl << "NFFT_plan::wipe: unable to get device no";
    return false;
  }

  if( device != old_device && cudaSetDevice(device) != cudaSuccess) {
    cerr << endl << "NFFT_plan::wipe: unable to set device no";
    return false;
  }

  if( mode==NFFT_WIPE_ALL && initialized ){
    deapodization_filter.reset();
    initialized = false;
  }
    
  if( preprocessed_NFFT_H ){
    delete tuples_last;
    delete bucket_begin;
    delete bucket_end;
  }

  if( preprocessed_NFFT || preprocessed_NFFT_H ){
    delete trajectory_positions;
    preprocessed_NFFT = preprocessed_NFFT_H = false;
  }

  if( device != old_device && cudaSetDevice(old_device) != cudaSuccess) {
    cerr << endl << "NFFT_plan::wipe: unable to restore device no";
    return false;
  }

  return true;
}

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::compute_beta()
{	
  // Compute Kaiser-Bessel beta paramter according to the formula provided in 
  // Beatty et. al. IEEE TMI 2005;24(6):799-808.

  beta = (M_PI*std::sqrt((W*W)/(alpha*alpha)*(alpha-0.5)*(alpha-0.5)-0.8)); 

  return true;
}

//
// Grid fictitious trajectory with a single sample at the origin
//

template<class REAL, unsigned int D> __global__ void
compute_deapodization_filter_kernel( typename uintd<D>::Type matrix_size_os, typename reald<REAL,D>::Type matrix_size_os_real, 
				     REAL W, REAL half_W, REAL one_over_W, REAL beta, typename complext<REAL>::Type*image_os )
{

  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int num_elements = prod(matrix_size_os);

  if( idx <num_elements ){

    // Compute weight from Kaiser-Bessel filter
    const typename uintd<D>::Type cell_pos = idx_to_co<D>(idx, matrix_size_os);

    // Sample position ("origin")
    const vector_td<REAL,D> sample_pos = get_half<REAL>()*matrix_size_os_real;

    // Calculate the distance between the cell and the sample
    vector_td<REAL,D> cell_pos_real = to_reald<REAL,unsigned int,D>(cell_pos);
    const typename reald<REAL,D>::Type delta = abs(sample_pos-cell_pos_real);

    // Compute convolution weight. 
    REAL weight; 
    REAL zero = get_zero<REAL>();
    vector_td<REAL,D> half_W_vec = to_vector_td<REAL,D>( half_W );

    if( weak_greater(delta, half_W_vec ) )
      weight = zero;
    else{ 
      weight = KaiserBessel<REAL>( delta, matrix_size_os_real, one_over_W, beta );
      //if( !isfinite(weight) )
      //weight = zero;
    }
    
    // Output weight
   typename complext<REAL>::Type result;
    result.vec[0] = weight; 
    result.vec[1] = zero;
    image_os[idx] = result;
  }
}

//
// Function to calculate the deapodization filter
//

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::compute_deapodization_filter()
{
  std::vector<unsigned int> tmp_vec_os = uintd_to_vector<D>(matrix_size_os);
  deapodization_filter = boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> >
    (cuNDArray<typename complext<REAL>::Type>::allocate(&tmp_vec_os));
  
  vector_td<REAL,D> matrix_size_os_real = to_reald<REAL,unsigned int, D>(matrix_size_os);
  
  // Be optimistic
  bool success = true;

  // Find dimensions of grid/blocks.

  dim3 dimBlock( 256 );
  dim3 dimGrid( (prod(matrix_size_os)+dimBlock.x-1)/dimBlock.x );

  // Invoke kernel
  compute_deapodization_filter_kernel<REAL,D><<<dimGrid, dimBlock>>> 
    ( matrix_size_os, matrix_size_os_real, W, get_half<REAL>()*W, reciprocal<REAL>(W), beta, deapodization_filter->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();
  
  // FFT
  if( success )
    success = fft( deapodization_filter.get(), NFFT_BACKWARDS, false );
  
  // Reciprocal
  if( success )
    cuNDA_reciprocal<typename complext<REAL>::Type>( deapodization_filter.get() );

  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::compute_NFFT( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image )
{
  // private method - no consistency check. We trust in ourselves.

  bool success = true;

  // Convolution
  if( success )
    success = deapodize( image );
    
  // FFT
  if( success )
    success = fft( image, NFFT_FORWARDS ); 

  // Deapodization
  if( success )
    success = convolve( samples, image, 0x0, NFFT_FORWARDS );
  
  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::compute_NFFT_H( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image )
{
  // private method - no consistency check. We trust in ourselves.

  bool success = true;

  // Convolution
  if( success )
    success = convolve( samples, image, 0x0, NFFT_BACKWARDS );
  
  // FFT
  if( success )
    success = fft( image, NFFT_BACKWARDS );
  
  // Deapodization
  if( success )
    success = deapodize( image );
  
  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::convolve_NFFT( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, bool accumulate )
{
  // private method - no consistency check. We trust in ourselves.

  if( !initialize_static_variables() ){
    cout << endl << "NFFT_plan::convolve_NFFT: unable to query device properties." << endl ;
    return false;
  }

  // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
  if( !((warp_size[device] & (warp_size[device]-1)) == 0 ) ){
    printf("\nError: on unsupported hardware (warpSize is not a power of two).\n");
    return false;
  }
  
  unsigned int num_batches = 1;
  for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
    num_batches *= image->get_size(d);
  num_batches /= number_of_frames;

  /*
    Setup grid and threads
  */

  size_t threads_per_block;
  unsigned int max_coils;

  switch(D){
    
  case 2:
    if( major[device] == 1 ){
      threads_per_block = NFFT_THREADS_PER_2D_KERNEL_1x;
      max_coils = NFFT_MAX_COILS_2D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_THREADS_PER_2D_KERNEL_2x;
      max_coils = NFFT_MAX_COILS_2D_COMPUTE_2x;
    }
    break;
    
  case 3:
    if( major[device] == 1 ){
      threads_per_block = NFFT_THREADS_PER_3D_KERNEL_1x;
      max_coils = NFFT_MAX_COILS_3D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_THREADS_PER_3D_KERNEL_2x;
      max_coils = NFFT_MAX_COILS_3D_COMPUTE_2x;
    }
    break;
    
  case 4:
    if( major[device] == 1 ){
      threads_per_block = NFFT_THREADS_PER_4D_KERNEL_1x;
      max_coils = NFFT_MAX_COILS_4D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_THREADS_PER_4D_KERNEL_2x;
      max_coils = NFFT_MAX_COILS_4D_COMPUTE_2x;
    }
    break;
    
  default:
    printf("\nNFFT Convolution: Dimensionality not supported!\n");
    return false;
  }
  
  // We can (only) convolve max_coils batches per run due to shared memory issues. 
  unsigned int domain_size_coils_desired = num_batches;
  unsigned int num_repetitions = domain_size_coils_desired/(max_coils+1) + 1;
  unsigned int domain_size_coils = (num_repetitions==1) ? domain_size_coils_desired : max_coils;
  unsigned int domain_size_coils_tail = (num_repetitions==1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions-1)*domain_size_coils;

  // Block and Grid dimensions
  dim3 dimBlock( (unsigned int)threads_per_block ); 
  dim3 dimGrid( (number_of_samples+dimBlock.x-1)/dimBlock.x, number_of_frames );

  // Calculate how much shared memory to use per thread
  size_t bytes_per_thread = domain_size_coils * sizeof( vector_td<REAL,D> );
  size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof( vector_td<REAL,D> );

  /*
    Invoke kernel
  */

  unsigned int double_warp_size_power=0;
  unsigned int __tmp = warp_size[device]<<1;
  while(__tmp!=1){
    __tmp>>=1;
    double_warp_size_power++;
  }
  
  vector_td<REAL,D> matrix_size_os_real = to_reald<REAL,unsigned int,D>( matrix_size_os );

  for( unsigned int repetition = 0; repetition<num_repetitions; repetition++ ){
    NFFT_convolve_kernel<REAL,D>
      <<<dimGrid, dimBlock, (repetition==num_repetitions-1) ? dimBlock.x*bytes_per_thread_tail : dimBlock.x*bytes_per_thread>>>
      ( alpha, beta, W, matrix_size_os, matrix_size_wrap, number_of_samples, 
	(repetition==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
	raw_pointer_cast(&(*trajectory_positions)[0]), 
	image->get_data_ptr()+repetition*prod(matrix_size_os)*number_of_frames*domain_size_coils,
	samples->get_data_ptr()+repetition*number_of_samples*number_of_frames*domain_size_coils, 
	double_warp_size_power, get_half<REAL>()*W, reciprocal(W), accumulate, matrix_size_os_real );
  }
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::convolve_NFFT_H( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, bool accumulate )
{

  // private method - no consistency check. We trust in ourselves.

  if( !initialize_static_variables() ){
    cout << endl << "NFFT_plan::convolve_NFFT: unable to query device properties." << endl ;
    return false;
  }

  // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
  if( !((warp_size[device] & (warp_size[device]-1)) == 0 ) ){
    printf("\nError: on unsupported hardware (warpSize is not a power of two).\n");
    return false;
  }
  
  unsigned int num_batches = 1;
  for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
    num_batches *= image->get_size(d);
  num_batches /= number_of_frames;

  /*
    Setup grid and threads
  */

  size_t threads_per_block;
  unsigned int max_coils;

  switch(D){

  case 2:
    if( major[device] == 1 ){
      threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_1x;
      max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_2x;
      max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_2x;
    }
    break;
      
  case 3:
    if( major[device] == 1 ){
      threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_1x;
      max_coils = NFFT_H_MAX_COILS_3D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_2x;
      max_coils = NFFT_H_MAX_COILS_3D_COMPUTE_2x;
    }
    break;
    
  case 4:
    if( major[device] == 1 ){
      threads_per_block = NFFT_H_THREADS_PER_4D_KERNEL_1x;
      max_coils = NFFT_H_MAX_COILS_4D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_H_THREADS_PER_4D_KERNEL_2x;
      max_coils = NFFT_H_MAX_COILS_4D_COMPUTE_2x;
    }
    break;
    
  default:
    printf("\nNFFT Convolution: Dimensionality not supported!\n");
    return false;
  }
  
  // We can (only) convolve domain_size_coils batches per run due to shared memory issues. 
  unsigned int domain_size_coils_desired = num_batches;
  unsigned int num_repetitions = domain_size_coils_desired/(max_coils+1) + 1;
  unsigned int domain_size_coils = (num_repetitions==1) ? domain_size_coils_desired : max_coils;
  unsigned int domain_size_coils_tail = (num_repetitions==1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions-1)*domain_size_coils;
  
  // Block and Grid dimensions
  dim3 dimBlock( (unsigned int)threads_per_block ); 
  dim3 dimGrid((prod(matrix_size_os+matrix_size_wrap)+dimBlock.x-1)/dimBlock.x, number_of_frames );
  
  // Calculate how much shared memory to use per thread
  size_t bytes_per_thread = domain_size_coils * sizeof( vector_td<REAL,D> );
  size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof( vector_td<REAL,D> );

  vector<unsigned int> vec_dims = uintd_to_vector<D>(matrix_size_os+matrix_size_wrap); 
  if( number_of_frames > 1 )
    vec_dims.push_back(number_of_frames);
  if( num_batches > 1 ) 
    vec_dims.push_back(num_batches);
  
  // Allocate and clear temporary image that includes a wrapping zone
  cuNDArray<typename complext<REAL>::Type> *_tmp = new cuNDArray<typename complext<REAL>::Type>; _tmp->create(&vec_dims); 
  
  if( !_tmp ){
    cout << endl << "Memory allocation failed before convolution." << endl;
    return false;
  }
  
  unsigned int double_warp_size_power=0, __tmp = warp_size[device]<<1;
  while(__tmp!=1){
    __tmp>>=1;
    double_warp_size_power++;
  }

  vector_td<REAL,D> matrix_size_os_real = to_reald<REAL,unsigned int,D>( matrix_size_os );

  for( unsigned int i = 0; i<num_repetitions; i++ ){
    NFFT_H_convolve_kernel<REAL,D>
      <<<dimGrid, dimBlock, (i==num_repetitions-1) ? dimBlock.x*bytes_per_thread_tail : dimBlock.x*bytes_per_thread>>>
      ( alpha, beta, W, matrix_size_os+matrix_size_wrap, number_of_samples, 
	(i==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
	raw_pointer_cast(&(*trajectory_positions)[0]), 
	_tmp->get_data_ptr()+i*prod(matrix_size_os+matrix_size_wrap)*number_of_frames*domain_size_coils,
	samples->get_data_ptr()+i*number_of_samples*number_of_frames*domain_size_coils, 
	raw_pointer_cast(&(*tuples_last)[0]), raw_pointer_cast(&(*bucket_begin)[0]), raw_pointer_cast(&(*bucket_end)[0]),
	double_warp_size_power, get_half<REAL>()*W, reciprocal(W), matrix_size_os_real );
  }
  
  CHECK_FOR_CUDA_ERROR();
  
  bool success = image_wrap( _tmp, image, accumulate );
 
  delete _tmp;

  return success;
}


// Image wrap kernels

template<class REAL, unsigned int D> __global__ void
image_wrap_kernel( typename uintd<D>::Type matrix_size_os, typename uintd<D>::Type matrix_size_wrap, bool accumulate,
		  typename complext<REAL>::Type*in,typename complext<REAL>::Type*out )
{
  const unsigned int number_of_images = blockDim.z;
  const unsigned int num_elements_per_image_src = prod(matrix_size_os+matrix_size_wrap);

  unsigned int idx;

  if( number_of_images == 1 )
    idx = (blockIdx.x/(matrix_size_os.vec[0]/blockDim.x))*blockDim.y*matrix_size_os.vec[0] + 
      (blockIdx.x%(matrix_size_os.vec[0]/blockDim.x))*blockDim.x + threadIdx.y*matrix_size_os.vec[0] + threadIdx.x;
  else
    idx = blockIdx.x*blockDim.x + threadIdx.x;

  const unsigned int image_offset_src = (threadIdx.z*gridDim.y+blockIdx.y)*num_elements_per_image_src;
  
  const typename uintd<D>::Type co = idx_to_co<D>(idx, matrix_size_os);
  const typename uintd<D>::Type half_wrap = matrix_size_wrap>>1;
  
  // Make "boolean" vectors denoting whether wrapping needs to be performed in a given direction (forwards/backwards)
  typename uintd<D>::Type B_l = vector_less( co, half_wrap ); 
  typename uintd<D>::Type B_r = vector_greater_equal( co, matrix_size_os-half_wrap ); 
  
  typename complext<REAL>::Type result = in[co_to_idx<D>(co+half_wrap, matrix_size_os+matrix_size_wrap) + image_offset_src];

  if( sum(B_l+B_r) > 0 ){
    
    // Fold back the wrapping zone onto the image ("periodically")
    //
    // There is 2^D-1 ways to pick combinations of dimensions in D-dimensionsal space, e.g. 
    // 
    //  { x, y, xy } in 2D
    //  { x, y, x, xy, xz, yz, xyz } in 3D
    //
    // Every "letter" in each combination provides two possible wraps (eiher end of the dimension)
    // 
    // For every 2^D-1 combinations DO
    //   - find the number of dimensions, d, in the combination
    //   - create 2^(d) stride vectors and test for wrapping using the 'B'-vectors above.
    //   - accumulate the contributions
    // 
    //   The following code represents dimensions as bits in a char.
    //
    
    for( unsigned char combination = 1; combination < (1<<D); combination++ ){
    
      // Find d
      unsigned char d = 0;
      for( unsigned char i=0; i<D; i++ )
	d += ((combination & (1<<i)) > 0 );
       
      // Create stride vector for each wrapping test
      for( unsigned char s = 0; s < (1<<d); s++ ){
        
	// Target for stride
	typename intd<D>::Type stride;
	char wrap_requests = 0;
	char skipped_dims = 0;
	
	// Fill dimensions of the stride
	for( unsigned char i=1; i<D+1; i++ ){
    
	  // Is the stride dimension present in the current combination?
	  if( i & combination ){
    
	    // A zero bit in s indicates "check for left wrap" and a one bit is interpreted as "check for right wrap" 
	    // ("left/right" for the individual dimension meaning wrapping on either side of the dimension).
    
	    if( i & (s<<(skipped_dims)) ){
	      if( B_r.vec[i-1] ){ // Wrapping required 
		set( stride, i-1, -1 );
		wrap_requests++;
	      }
	      else
		set( stride, i-1, 0 );
	    }
	    else{ 
	      if( B_l.vec[i-1] ){ // Wrapping required 
		set( stride, i-1, +1 );
		wrap_requests++;
	      }
	      else
		set( stride, i-1, 0 );
	    }
	  }
	  else{
	    // Do not test for wrapping in dimension 'i-1' (for this combination)
	    set( stride, i-1, 0 );
	    skipped_dims++;
	  }
	}
	
	// Now it is time to do the actual wrapping (if needed)
	if( wrap_requests == d ){
	  typename intd<D>::Type src_co_int = to_intd(co+half_wrap);
	  typename intd<D>::Type matrix_size_os_int = to_intd(matrix_size_os);
	  typename intd<D>::Type co_offset_int = src_co_int + component_wise_mul<int,D>(stride,matrix_size_os_int);
	  typename uintd<D>::Type co_offset = to_uintd(co_offset_int);
	  result += in[co_to_idx<D>(co_offset, matrix_size_os+matrix_size_wrap) + image_offset_src];
	  break; // only one stride per combination can contribute (e.g. one edge, one corner)
	} 
      } 
    }
  }
  
  // Output
  const unsigned int image_offset_tgt = (threadIdx.z*gridDim.y+blockIdx.y)*prod(matrix_size_os);
  if( accumulate ) result += out[idx+image_offset_tgt];
  out[idx+image_offset_tgt] = result;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::image_wrap( cuNDArray<typename complext<REAL>::Type> *source, cuNDArray<typename complext<REAL>::Type> *target, bool accumulate )
{
  unsigned int num_batches = 1;
  for( unsigned int d=D; d<source->get_number_of_dimensions(); d++ )
    num_batches *= source->get_size(d);
  num_batches /= number_of_frames;

  // Set dimensions of grid/blocks.
  unsigned int bdim = 16;
  dim3 dimBlock( bdim, (num_batches==1) ? bdim : 1, num_batches );
  dim3 dimGrid( prod(matrix_size_os)/(dimBlock.x*dimBlock.y), number_of_frames );

  // Safety check
  if( (matrix_size_os.vec[0]%bdim) || (num_batches==1 && (matrix_size_os.vec[1]%bdim)) ) {
    cout << endl << "dimensions mismatch in wrapping setup." << endl;
    return false;
  }

  // Invoke kernel
  image_wrap_kernel<REAL,D><<<dimGrid, dimBlock>>>( matrix_size_os, matrix_size_wrap, accumulate, source->get_data_ptr(), target->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}	

//
// Template instantion
//

template class EXPORTGPUNFFT NFFT_plan< float, 2 >;
template class EXPORTGPUNFFT NFFT_plan< double, 2 >;

template class EXPORTGPUNFFT NFFT_plan< float, 3 >;
template class EXPORTGPUNFFT NFFT_plan< double, 3 >;

template class EXPORTGPUNFFT NFFT_plan< float, 4 >;
template class EXPORTGPUNFFT NFFT_plan< double, 4 >;
