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
#include "hoNDArray_fileio.h"


// Includes - our own code
#include "NFFT.hcu"
#include "cuNDFFT.h"

#include "cuNDArray.h"
#include "ndarray_device_utilities.hcu"

#include "uintd_operators.hcu"
#include "uintd_utilities.hcu"

#include "floatd_operators.hcu"
#include "floatd_utilities.hcu"

#include "check_CUDA.h"

// Includes - CUDA
#include <math_constants.h>
#include <cuComplex.h>
#include <cufft.h>
#include <cublas.h>

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

// TODO: make d x generation dimensional struct to contain MAX_COILS and THREADS_PER_KERNEL for the NFFT and NFFT_H respectively

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

#include "NFFT_kernel.cu"
#include "NFFT_H_kernel.cu"
#include "KaiserBessel_kernel.cu"
#include "NFFT_preprocess_kernel.cu"

//
// Public class methods
//

template<class UINTd, class REALd, class REAL, class NDTYPE> 
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::NFFT_plan()
{
  // Minimal initialization
  barebones();
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::NFFT_plan( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, REAL W, int device )
{
  // Minimal initialization
  barebones();

  // Setup plan
  setup( matrix_size, matrix_size_os, fixed_dims, W, device );
  
  if( !initialized )
    cout << endl << "Initialization of the plan failed." << endl;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::NFFT_plan( NFFT_plan<UINTd, REALd, REAL, NDTYPE> *plan )
{
  matrix_size = plan->matrix_size;
  matrix_size_os = plan->matrix_size_os;
  matrix_size_wrap = plan->matrix_size_wrap;
  fixed_dims = plan->fixed_dims;
  non_fixed_dims = plan->non_fixed_dims;
  d = plan->d;
  W = plan->W;
  alpha = plan->alpha;
  beta = plan->beta;
  number_of_samples = plan->number_of_samples;
  
  if (cudaGetDevice(&device) != cudaSuccess) {
    cerr << "NFFT_plan::NFFT_plan: unable to get device no" << endl;
  }  

  if(plan->initialized){
    *deapodization_filter = *(plan->deapodization_filter);
    initialized = true;
  }

  preprocessed_NFFT = plan->preprocessed_NFFT;
  preprocessed_NFFT_H = plan->preprocessed_NFFT_H;
   
  if( plan->preprocessed_NFFT || plan->preprocessed_NFFT_H ){
      *trajectory_positions = *(plan->trajectory_positions);
    }
 
  if( plan->preprocessed_NFFT_H ){
    *tuples_last = *(plan->tuples_last);
    *bucket_begin = *(plan->bucket_begin);
    *bucket_end = *(plan->bucket_end);
  }
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::~NFFT_plan()
{
  wipe();
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
bool NFFT_plan<UINTd, REALd, REAL, NDTYPE>::setup( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, REAL W, int _device )
{	
  // Free memory
  wipe();

  //
  // Check if the device is valid
  //

  if( _device<0 )
    cudaGetDevice( (int*)&device );
  else
    device = _device;

  cudaDeviceProp deviceProp;  
  cudaGetDeviceProperties( &deviceProp, device );

  unsigned int warp_size = deviceProp.warpSize;
  
  //
  // Check input against certain requirements
  //
  
  if( sum(matrix_size%warp_size) || sum(matrix_size_os%warp_size) ){
    cout << endl << "Illegal matrix size for the NFFT plan (not a multiple of the warp size)" << endl;
    return false;
  }

  //
  // Setup private variables
  //

  d = sizeof(UINTd)/sizeof(unsigned int);

  this->matrix_size = matrix_size;
  this->matrix_size_os = matrix_size_os;
  this->fixed_dims = fixed_dims;

  for( unsigned int i=0; i<d; i++ ){

    if( ((unsigned int*)(&fixed_dims))[i] >= 1 ){
      ((unsigned int*)(&non_fixed_dims))[i] = 0;
    }
    else{
      ((unsigned int*)(&non_fixed_dims))[i] = 1;
    }
  }

  REAL one; get_one(one);
  UINTd ones; uint_to_uintd( 1, 1, ones );  
  REAL W_half = half(W);
  REALd W_vec; real_to_reald( W_half, W_half, W_vec );

  matrix_size_wrap = hreald_to_uintd(ceil(W_vec))<<1; // TODO: check if padding to warp_size changes performance
  
  alpha = (REAL) matrix_size_os.x / (REAL) matrix_size.x;
  
  if( alpha < one ){
    cout << endl << "Illegal oversampling ratio suggested" << endl;
    return false;
  }

  this->W = W;

  unsigned int fracX = matrix_size_os.x / matrix_size.x;
  unsigned int moduX = matrix_size_os.x % matrix_size.x;

  for( unsigned int dim=0; dim<d; dim++){

    if( !(((unsigned int*)(&fixed_dims))[dim]) && ( W > ((unsigned int*)(&matrix_size))[dim])){
    }

    if( (((unsigned int*)(&fixed_dims))[dim]) > 1 )
      ((unsigned int*)(&fixed_dims))[dim] = 1; // enforce fixed_dims as "boolean style"

    if( (((unsigned int*)(&fixed_dims))[dim]) ){
      ((unsigned int*)(&matrix_size_wrap))[dim] = 0;
    }
    
    if( !(((unsigned int*)(&fixed_dims))[dim]) &&
	((((unsigned int*)(&matrix_size_os))[dim]/(((unsigned int*)(&matrix_size))[dim])) != fracX ||
	 (((unsigned int*)(&matrix_size_os))[dim]%(((unsigned int*)(&matrix_size))[dim])) != moduX) ){
      cout << endl << "Oversampling ratio is not constant between dimensions" << endl;
      return false;
    }
  }
  
  // Compute Kaiser-Bessel beta
  bool success = 
    compute_beta();
  
  //cout << endl << "beta: " << beta << endl;

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
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }
  
  return success;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
bool NFFT_plan<UINTd, REALd, REAL, NDTYPE>::preprocess( cuNDArray<REALd> *trajectory, bool forwards_only )
{
  if( trajectory->get_device() != device ){
    cout << endl << "NFFT_plan::preprocess: device mismatch." << endl;
    return false;
  }

  wipe( true );

  if( !trajectory || !trajectory->get_number_of_dimensions() == 1 ){
    cout << endl << "NFFT_plan preprocessing: illegal trajecotry ndarray." << endl;
    return false;
  }
  
  if( !initialized ){
    cout << endl << "NFFT_plan preprocessing: Unable to proceed without setup." << endl;
    return false;
  }

  number_of_samples = trajectory->get_size(0);

  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to get device no" << endl;
  }  
  if( device != device_no_old && cudaSetDevice(device) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }  
  
  // Make Thrust device vector of trajectory and samples
  device_vector<REALd> trajectory_positions_in( device_pointer_cast<REALd>(trajectory->get_data_ptr()), device_pointer_cast<REALd>(trajectory->get_data_ptr()+number_of_samples) );
  trajectory_positions = new device_vector<REALd>( number_of_samples );

  CHECK_FOR_CUDA_ERROR();

  // convert input trajectory in [-1/2;1/2] to [0;matrix_size_os]
  thrust::transform( trajectory_positions_in.begin(), trajectory_positions_in.end(), trajectory_positions->begin(), trajectory_scale<REALd>(huintd_to_reald(matrix_size_os), huintd_to_reald((matrix_size_os+matrix_size_wrap)>>1)) );

  if( !forwards_only ){
    
    // allocate storage for and compute temporary prefix-sum variable (#cells influenced per sample)
    device_vector<unsigned int> c_p_s(number_of_samples);
    device_vector<unsigned int> c_p_s_ps(number_of_samples);
    
    CHECK_FOR_CUDA_ERROR();
    
    REAL half_W = half(W);
    thrust::plus<unsigned int> binary_op;
    thrust::transform(trajectory_positions->begin(), trajectory_positions->end(), c_p_s.begin(), compute_num_cells_per_sample<REALd, REAL>(d,half_W));
    inclusive_scan( c_p_s.begin(), c_p_s.end(), c_p_s_ps.begin(), binary_op ); // prefix sum
    
    // Build the vector of (grid_idx, sample_idx) tuples. Actually kept in two seperate vectors.
    unsigned int num_pairs = c_p_s_ps.back();
    thrust::device_vector<unsigned int> tuples_first = device_vector<unsigned int>(num_pairs);
    tuples_last = new device_vector<unsigned int>(num_pairs);
    
    CHECK_FOR_CUDA_ERROR();
    
    // Fill tuple vector
    write_pairs<UINTd, REALd, REAL>( matrix_size_os, matrix_size_wrap, number_of_samples, W, 
				     raw_pointer_cast(&(*trajectory_positions)[0]), raw_pointer_cast(&c_p_s_ps[0]), 
				     raw_pointer_cast(&tuples_first[0]), raw_pointer_cast(&(*tuples_last)[0]) );
    
    // Sort by grid indices
    sort_by_key(tuples_first.begin(), tuples_first.end(), tuples_last->begin() );
    
    // each bucket_begin[i] indexes the first element of bucket i's list of points
    // each bucket_end[i] indexes one past the last element of bucket i's list of points
    bucket_begin = new device_vector<unsigned int>(prod(matrix_size_os+matrix_size_wrap));
    bucket_end = new device_vector<unsigned int>(prod(matrix_size_os+matrix_size_wrap));
    
    CHECK_FOR_CUDA_ERROR();
    
    // find the beginning of each bucket's list of points
    counting_iterator<unsigned int> search_begin(0);
    lower_bound(tuples_first.begin(), tuples_first.end(), search_begin,	search_begin + prod(matrix_size_os+matrix_size_wrap), bucket_begin->begin() );
    
    // find the end of each bucket's list of points
    upper_bound(tuples_first.begin(), tuples_first.end(), search_begin, search_begin + prod(matrix_size_os+matrix_size_wrap), bucket_end->begin() );
  }
  
  if( device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }
  
  preprocessed_NFFT = true;

  if( !forwards_only )
    preprocessed_NFFT_H = true;

  return true;
}
  

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::compute( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image, cuNDArray<REAL> *weights, NFFT_mode mode )
{  
  if( samples->get_device() != device || image->get_device() != device || (weights && weights->get_device() != device) ){
    cout << endl << "NFFT_plan::compute: device mismatch." << endl;
    return false;
  }

  // Validity checks
  
  unsigned char components;

  if( mode == NFFT_FORWARDS ) 
    components = NFFT_CONVOLUTION + NFFT_FFT + NFFT_DEAPODIZATION;

  else // backwards
    components = NFFT_H_CONVOLUTION + NFFT_FFT + NFFT_DEAPODIZATION;

  if( !check_consistency( samples, image, weights, components ) )
    return false;
  
  bool success;  
  cuNDArray<NDTYPE> *working_image = 0x0, *working_samples = 0x0;
  
  UINTd image_dims; 
  if( !cuNDA_fromVec( image->get_dimensions(), image_dims ) ){
    cout << "NFFT_plan::compute: image dimension undermatch the plan" << endl;
    return false;
  }
  
  bool oversampled_image = (image_dims==matrix_size_os); 
  unsigned int num_batches = (image->get_number_of_dimensions()>d) ? image->get_size(image->get_number_of_dimensions()-1) : 1;
  vector<unsigned int> vec_dims = cuNDA_toVec(matrix_size_os); 
  if( num_batches > 1 ) 
    vec_dims.push_back(num_batches);

  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to get device no" << endl;
  }  
  if( device != device_no_old && cudaSetDevice(device) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }  

  switch(mode){

  case NFFT_FORWARDS:
    
    if( !oversampled_image ){
      working_image = cuNDArray<NDTYPE>::allocate(vec_dims);
      cuNDA_expand_with_zero_fill<UINTd, REALd>( image, working_image );
    }
    else{
      working_image = image;
    }
    
    success = compute_NFFT( samples, working_image );

    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    break;
    
  case NFFT_BACKWARDS:

    // Density compensation
    if( weights ){
      working_samples = new cuNDArray<NDTYPE>(*samples);
      cuNDA_scale( weights, working_samples );
    }
    else{
      working_samples = samples;
    }
    
    if( !oversampled_image ){
      working_image = cuNDArray<NDTYPE>::allocate(vec_dims);
    }
    else{
      working_image = image;
    }

    success = compute_NFFT_H( working_samples, working_image );

    if( success ){
      if( !oversampled_image ){
	cuNDA_crop( (matrix_size_os-matrix_size)>>1, working_image, image );
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

  if( device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }

  return success;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::compute_iteration( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image, cuNDArray<REAL> *weights, NFFT_mode mode )
{  
  if( samples->get_device() != device || image->get_device() != device || (weights && weights->get_device() != device) ){
    cout << endl << "NFFT_plan::compute_iteration: device mismatch." << endl;
    return false;
  }

  // Validity checks
  
  unsigned char components = NFFT_CONVOLUTION + NFFT_H_CONVOLUTION + NFFT_FFT + NFFT_DEAPODIZATION;
  
  if( !check_consistency( samples, image, weights, components ) )
    return false;
  
  bool success;  
  cuNDArray<NDTYPE> *working_image = 0x0, *working_samples = 0x0;

  UINTd image_dims; 
  if( !cuNDA_fromVec( image->get_dimensions(), image_dims ) ){
    cout << "NFFT_plan::compute_iteration: image dimension undermatch the plan" << endl;
    return false;
  }
  
  bool oversampled_image = (image_dims==matrix_size_os); 
  unsigned int num_batches = (image->get_number_of_dimensions()>d) ? image->get_size(image->get_number_of_dimensions()-1) : 1;
  vector<unsigned int> vec_dims = cuNDA_toVec(matrix_size_os); 
  if( num_batches > 1 ) 
    vec_dims.push_back(num_batches);

  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to get device no" << endl;
  }  
  if( device != device_no_old && cudaSetDevice(device) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }  

  switch(mode){

  case NFFT_FORWARDS: // iteration from image

    if( !oversampled_image ){
      working_image = cuNDArray<NDTYPE>::allocate(vec_dims);
      cuNDA_expand_with_zero_fill<UINTd, REALd>( image, working_image );
    }
    else{
      working_image = image;
    }
    
    working_samples = samples;

    success = compute_NFFT( working_samples, working_image );
    
    if( success ){
      
      // Density compensation
      if( weights ){
	cuNDA_scale( weights, working_samples );
      }
      
      success = compute_NFFT_H( working_samples, working_image ); 
    }      
    
    if( success ){
      if( !oversampled_image ){
	cuNDA_crop( (matrix_size_os-matrix_size)>>1, working_image, image );
	delete working_image; working_image = 0x0;
      }
    }
    
    break;

  case NFFT_BACKWARDS: // iteration from samples

    // Density compensation
    if( weights ){
      working_samples = new cuNDArray<NDTYPE>(*samples);
      cuNDA_scale( weights, working_samples );
    }
    else{
      working_samples = samples;
    }
    
    if( !oversampled_image )
      working_image = cuNDArray<NDTYPE>::allocate(vec_dims);
    else
      working_image = image;

    success = compute_NFFT_H( working_samples, working_image );
    
    if( success ){
      
      cuNDA_zero_fill_border( matrix_size, working_image );
      
      success = compute_NFFT( samples, working_image );
    }
  
    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    if( weights ){
      delete working_samples; working_samples = 0x0;
    }
    
    break;
  };

  if( device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }

  return success;
}

template <class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::convolve( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image, cuNDArray<REAL> *weights, NFFT_mode mode )
{
  if( samples->get_device() != device || image->get_device() != device || (weights && weights->get_device() != device) ){
    cout << endl << "NFFT_plan::convolve: device mismatch." << endl;
    return false;
  }

  unsigned char components;

  if( mode == NFFT_FORWARDS ) 
    components = NFFT_CONVOLUTION;
  else
    components = NFFT_H_CONVOLUTION;
  
  if( !check_consistency( samples, image, weights, components ) )
    return false;
  
  bool success;  
  cuNDArray<NDTYPE> *working_image = 0x0, *working_samples = 0x0;
  
  UINTd image_dims; 
  if( !cuNDA_fromVec( image->get_dimensions(), image_dims ) ){
    cout << "NFFT_plan::convolve: image dimension undermatch the plan" << endl;
    return false;
  }
  
  bool oversampled_image = (image_dims==matrix_size_os); 
  unsigned int num_batches = (image->get_number_of_dimensions()>d) ? image->get_size(image->get_number_of_dimensions()-1) : 1;
  vector<unsigned int> vec_dims = cuNDA_toVec(matrix_size_os); 
  if( num_batches > 1 ) 
    vec_dims.push_back(num_batches);

  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to get device no" << endl;
  }  
  if( device != device_no_old && cudaSetDevice(device) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }  

  switch(mode){

  case NFFT_FORWARDS:
    
    if( !oversampled_image ){
      working_image = cuNDArray<NDTYPE>::allocate(vec_dims);
      cuNDA_expand_with_zero_fill<UINTd, REALd>( image, working_image );
    }
    else{
      working_image = image;
    }
    
    success = convolve_NFFT( samples, working_image );

    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    break;
    
  case NFFT_BACKWARDS:

    // Density compensation
    if( weights ){
      working_samples = new cuNDArray<NDTYPE>(*samples);
      cuNDA_scale( weights, working_samples );
    }
    else{
      working_samples = samples;
    }
    
    if( !oversampled_image ){
      working_image = cuNDArray<NDTYPE>::allocate(vec_dims);
    }
    else{
      working_image = image;
    }
    
    success = convolve_NFFT_H( working_samples, working_image );

    if( success ){
      if( !oversampled_image ){
	cuNDA_crop( (matrix_size_os-matrix_size)>>1, working_image, image );
      }
    }
    
    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    if( weights ){
      delete working_samples; working_samples = 0x0;
    }
    
    break;
  }

  if( device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }

  return success;
}


template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::FFT(cuNDArray<NDTYPE> *data, NFFT_mode mode, bool do_scale )
{
  if( data->get_device() != device ){
    cout << endl << "NFFT_plan::FFT: device mismatch." << endl;
    return false;
  }

  UINTd _dims_to_transform; counting_vec( _dims_to_transform );
  vector<unsigned int> dims_to_transform = cuNDA_toVec( _dims_to_transform );

  int res;
  if( mode == NFFT_FORWARDS ){
    res = cuNDFFT().fft( data, dims_to_transform );
  }
  else{
    res = cuNDFFT().ifft( data, dims_to_transform, do_scale );
  }

  if( res == 0 )
    return true;
  else
    return false;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::deapodize( cuNDArray<NDTYPE> *image )
{
  if( image->get_device() != device ){
    cout << endl << "NFFT_plan::FFT: device mismatch." << endl;
    return false;
  }

  unsigned char components;
  
  components = NFFT_FFT;
  
  if( !check_consistency( 0x0, image, 0x0, components ) )
    return false;
  
  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to get device no" << endl;
  }  
  if( device != device_no_old && cudaSetDevice(device) != cudaSuccess) {
    cerr << "NFFT_plan::setup: unable to set device" << endl;
  }  
  
 cuNDA_scale( deapodization_filter, image );
  
 if( device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
   cerr << "NFFT_plan::setup: unable to set device" << endl;
 }

  return true;
}

//
// Private class methods
//

template<class UINTd, class REALd, class REAL, class NDTYPE> bool 
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::check_consistency( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image, cuNDArray<REAL> *weights, unsigned char components )
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

  if( !((image->get_number_of_dimensions() == d) || (image->get_number_of_dimensions() == (d+1))) ){
    cout << endl << "NFFT_plan: Number of image dimensions mismatch the plan." << endl;
    return false;
  }    

  UINTd image_dims; cuNDA_fromVec( image->get_dimensions(), image_dims );
  bool oversampled_image = (image_dims==matrix_size_os);
  
  if( !((oversampled_image) ? (image_dims == matrix_size_os) : (image_dims == matrix_size) )){
    cout << endl << "NFFT_plan: Image dimensions mismatch." << endl;
    return false;
  }
  
  if( (components & NFFT_CONVOLUTION ) || (components & NFFT_H_CONVOLUTION ) ){
    
    if( !(samples->get_number_of_dimensions()==1 || samples->get_number_of_dimensions()==2) ){
      cout << endl << "NFFT_plan: Number of sample dimensions should be one or two: samples_per_recon (x batch_size)." << endl;
      return false;
    }
    
    if( !(samples->get_size(0)==number_of_samples) ){
      cout << endl << "NFFT_plan: Number of samples mismatch the preprocessed trajectory." << endl;
      return false;
    }
    
    if( samples->get_number_of_dimensions()==2 || image->get_number_of_dimensions() == d+1 ){
      if( !(samples->get_number_of_dimensions()==2) || !(image->get_number_of_dimensions() == d+1) || !(samples->get_size(1)==image->get_size(d)) ){
	cout << endl << "NFFT_plan: Number of batches mismatch between samples and image." << endl;
	return false;
      }
    }
  }
  
  if( components & NFFT_H_CONVOLUTION ){
    if( weights ){
 
      if( !(weights->get_number_of_dimensions()==1 || weights->get_number_of_dimensions()==2) ){
	cout << endl << "NFFT_plan: Number of weight dimensions should be one or two: samples_per_recon (x batch_size)." << endl;
	return false;
      }
      
      if( !(weights->get_size(0)==number_of_samples) ){
	cout << endl << "NFFT_plan: Number of weights mismatch the preprocessed trajectory." << endl;
	return false;
      }
    }
  }

  return true;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
bool NFFT_plan<UINTd, REALd, REAL, NDTYPE>::barebones()
{	
  // These are the fundamental booleans checked before accessing the various member pointers
  initialized = preprocessed_NFFT = preprocessed_NFFT_H = false;
  
  deapodization_filter = 0x0;

  // and specify the device
  if (cudaGetDevice(&device) != cudaSuccess) {
    cerr << "NFFT_plan::barebones:: unable to get device no" << endl;
    return false;
  }

  return true;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
void NFFT_plan<UINTd, REALd, REAL, NDTYPE>::wipe( bool preprocess_only )
{
  if( !preprocess_only && initialized ){
    delete deapodization_filter;
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
}

template<class UINTd, class REALd, class REAL, class NDTYPE> 
bool NFFT_plan<UINTd, REALd, REAL, NDTYPE>::compute_beta()
{	
  // Compute Kaiser-Bessel beta paramter according to the formula provided in 
  // Beatty et. al. IEEE TMI 2005;24(6):799-808.

  beta = (M_PI*sqrt((W*W)/(alpha*alpha)*(alpha-0.5)*(alpha-0.5)-0.8)); 

  return true;
}

//
// Grid fictitious trajectory with a single sample at the origin
//

template <class UINTd, class REALd, class REAL, class NDTYPE> __global__ void
compute_deapodization_filter_kernel( UINTd matrix_size_os, REALd matrix_size_os_real, UINTd fixed_dims, 
				     REAL W, REAL half_W, REAL one_over_W, REAL beta, NDTYPE *image_os )
{

  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int num_elements = prod(matrix_size_os);

  if( idx <num_elements ){

    // Compute weight from Kaiser-Bessel filter
    const UINTd cell_pos = idx_to_co(idx, matrix_size_os);

    // Sample position ("origin")
    const REALd sample_pos = half(matrix_size_os_real);

    // Calculate the distance between the cell and the sample
    const REALd delta = abs(sample_pos-uintd_to_reald(cell_pos));

    // Compute convolution weight. 
    REAL weight, zero;
    get_zero(zero);
    REALd half_W_vec; real_to_reald(half_W, half_W, half_W_vec );

    if( weak_greater(delta, half_W_vec ) )
      weight = zero;
    else{ 
      weight = KaiserBessel( delta, matrix_size_os_real, one_over_W, beta, fixed_dims );
      if( !isfinite(weight) )
	weight = zero;
    }
    
    // Output weight
    image_os[idx] =  make_realComplex( weight, zero );
  }
}

//
// Function to calculate the deapodization filter
//

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::compute_deapodization_filter()
{

  if( initialized && deapodization_filter ) 
    delete( deapodization_filter );
  
  deapodization_filter = cuNDArray<NDTYPE>::allocate(cuNDA_toVec(matrix_size_os));

  if( !deapodization_filter ){
    cout << endl << "cuNDArray allocation failed" << endl;
    return false;
  }
  
  // Be optimistic
  bool success = true;

  // Find dimensions of grid/blocks.

  dim3 dimBlock( 256 );
  dim3 dimGrid( (unsigned int) ceil((double)prod(matrix_size_os)/(double)dimBlock.x) );

  // Invoke kernel
  compute_deapodization_filter_kernel<UINTd, REALd, REAL, NDTYPE><<<dimGrid, dimBlock>>> 
    ( matrix_size_os, huintd_to_reald(matrix_size_os), fixed_dims, W, half(W), reciprocal(W), beta, deapodization_filter->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();
  
  // FFT
  if( success )
    success = FFT( deapodization_filter, NFFT_BACKWARDS );
  
  // Reciprocal
  if( success )
    cuNDA_reciprocal( deapodization_filter );

  if( !success ){
    delete deapodization_filter;
  }
  
  return success;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::compute_NFFT( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image )
{
  // private method - no consistency check. We trust in ourselves.

  bool success = true;

  // Convolution
  if( success )
    success = deapodize( image );
  
  // FFT
  if( success )
    success = FFT( image, NFFT_FORWARDS );
  
  // Deapodization
  if( success )
    success = convolve( samples, image, 0x0, NFFT_FORWARDS );
  
  return success;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::compute_NFFT_H( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image )
{
  // private method - no consistency check. We trust in ourselves.

  bool success = true;

  // Convolution
  if( success )
    success = convolve( samples, image, 0x0, NFFT_BACKWARDS );
  
  // FFT
  if( success )
    success = FFT( image, NFFT_BACKWARDS );
  
  // Deapodization
  if( success )
    success = deapodize( image );
  
  return success;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::convolve_NFFT( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image )
{
  // private method - no consistency check. We trust in ourselves.

  // Get device properties
  cudaDeviceProp deviceProp;  
  cudaGetDeviceProperties( &deviceProp, device );

  const unsigned int warp_size = deviceProp.warpSize;
  const unsigned int generation = deviceProp.major;

  // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
  if( !((warp_size & (warp_size-1)) == 0 ) ){
    printf("\nError: on unsupported hardware (warpSize is not a power of two).\n");
    return false;
  }
  
  /*
    Setup grid and threads
  */

  size_t threads_per_block;
  unsigned int max_coils;

  switch( sizeof(UINTd) )
    {

    case sizeof(uint2):
      if( generation == 1 ){
	threads_per_block = NFFT_THREADS_PER_2D_KERNEL_1x;
	max_coils = NFFT_MAX_COILS_2D_COMPUTE_1x;
      }
      else{
	threads_per_block = NFFT_THREADS_PER_2D_KERNEL_2x;
	max_coils = NFFT_MAX_COILS_2D_COMPUTE_2x;
      }
      break;
      
    case sizeof(uint3):
      if( generation == 1 ){
	threads_per_block = NFFT_THREADS_PER_2D_KERNEL_1x;
	max_coils = NFFT_MAX_COILS_2D_COMPUTE_1x;
      }
      else{
	threads_per_block = NFFT_THREADS_PER_2D_KERNEL_2x;
	max_coils = NFFT_MAX_COILS_2D_COMPUTE_2x;
      }
      break;

    case sizeof(uint4):
      if( generation == 1 ){
	threads_per_block = NFFT_THREADS_PER_3D_KERNEL_1x;
	max_coils = NFFT_MAX_COILS_3D_COMPUTE_1x;
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

  // We can (only) convolve domain_size_coils batches per run due to shared memory issues. 
  unsigned int domain_size_coils_desired = (samples->get_number_of_dimensions()==1) ? 1 : samples->get_size(1);
  unsigned int num_repetitions = domain_size_coils_desired/(max_coils+1) + 1;
  unsigned int domain_size_coils = (num_repetitions==1) ? domain_size_coils_desired : max_coils;
  unsigned int domain_size_coils_tail = (num_repetitions==1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions-1)*domain_size_coils;

  // Grid dimensions	
  unsigned int gridDimX = (unsigned int) ceil((double)number_of_samples/(double)threads_per_block);

  // Block and Grid dimensions
  dim3 dimBlock( (unsigned int)threads_per_block ); 
  dim3 dimGrid( gridDimX );

  // Calculate how much shared memory to use per thread
  size_t bytes_per_thread = domain_size_coils * sizeof(NDTYPE);
  size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(NDTYPE);

  /*
    Invoke kernel
  */

  unsigned int double_warp_size_power=0;
  unsigned int __tmp = warp_size<<1;
  while(__tmp!=1){
    __tmp>>=1;
    double_warp_size_power++;
  }
  
  for( unsigned int batch = 0; batch<num_repetitions; batch++ ){
    NFFT_convolve_kernel<UINTd,REALd,REAL,NDTYPE>
      <<<dimGrid, dimBlock, (batch==num_repetitions-1) ? prod(dimBlock)*bytes_per_thread_tail : prod(dimBlock)*bytes_per_thread>>>
      ( alpha, beta, W, matrix_size_os, matrix_size_wrap, fixed_dims, number_of_samples, 
	(batch==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
	raw_pointer_cast(&(*trajectory_positions)[0]), 
	&(image->get_data_ptr()[batch*prod(matrix_size_os)*domain_size_coils]), &(samples->get_data_ptr()[batch*number_of_samples*domain_size_coils]), 
	double_warp_size_power, half(W), reciprocal(W), huintd_to_reald(matrix_size_os), non_fixed_dims );
  }
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::convolve_NFFT_H( cuNDArray<NDTYPE> *samples, cuNDArray<NDTYPE> *image )
{

  // private method - no consistency check. We trust in ourselves.

  // Get device properties
  cudaDeviceProp deviceProp;  
  cudaGetDeviceProperties( &deviceProp, device );

  const unsigned int warp_size = deviceProp.warpSize;
  const unsigned int generation = deviceProp.major;

  // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
  if( !((warp_size & (warp_size-1)) == 0 ) ){
    printf("\nError: on unsupported hardware (warpSize is not a power of two).\n");
    return false;
  }
  
  /*
    Setup grid and threads
  */

  size_t threads_per_block;
  unsigned int max_coils;

  switch( sizeof(UINTd) )
    {

    case sizeof(uint2):
      if( generation == 1 ){
	threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_1x;
	max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_1x;
      }
      else{
	threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_2x;
	max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_2x;
      }
      break;
      
    case sizeof(uint3):
      if( generation == 1 ){
	threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_1x;
	max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_1x;
      }
      else{
	threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_2x;
	max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_2x;
      }
      break;

    case sizeof(uint4):
      if( generation == 1 ){
	threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_1x;
	max_coils = NFFT_H_MAX_COILS_3D_COMPUTE_1x;
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
  unsigned int domain_size_coils_desired = (samples->get_number_of_dimensions()==1) ? 1 : samples->get_size(1);
  unsigned int num_repetitions = domain_size_coils_desired/(max_coils+1) + 1;
  unsigned int domain_size_coils = (num_repetitions==1) ? domain_size_coils_desired : max_coils;
  unsigned int domain_size_coils_tail = (num_repetitions==1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions-1)*domain_size_coils;
  
  // Grid dimensions	
  unsigned int gridDimX = (unsigned int) ceil((double)prod(matrix_size_os+matrix_size_wrap)/(double)threads_per_block);

  // Block and Grid dimensions
  dim3 dimBlock( (unsigned int)threads_per_block ); 
  dim3 dimGrid( gridDimX );

  // Calculate how much shared memory to use per thread
  size_t bytes_per_thread = domain_size_coils * sizeof(NDTYPE);
  size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof(NDTYPE);

  vector<unsigned int> vec_dims = cuNDA_toVec(matrix_size_os+matrix_size_wrap); 
  if( domain_size_coils_desired > 1 ) 
    vec_dims.push_back(domain_size_coils_desired);
  
  // Allocate and clear temporary image that includes a wrapping zone
  cuNDArray<NDTYPE> *_tmp = cuNDArray<NDTYPE>::allocate(vec_dims);
  
  if( !_tmp ){
    cout << endl << "Memory allocation failed before convolution." << endl;
    return false;
  }
  
  unsigned int double_warp_size_power=0, __tmp = warp_size<<1;
  while(__tmp!=1){
    __tmp>>=1;
    double_warp_size_power++;
  }
  
  for( unsigned int batch = 0; batch<num_repetitions; batch++ ){
    NFFT_H_convolve_kernel<UINTd,REALd,REAL,NDTYPE>
      <<<dimGrid, dimBlock, (batch==num_repetitions-1) ? prod(dimBlock)*bytes_per_thread_tail : prod(dimBlock)*bytes_per_thread>>>
      ( alpha, beta, W, matrix_size_os+matrix_size_wrap, fixed_dims, number_of_samples, 
	(batch==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
	raw_pointer_cast(&(*trajectory_positions)[0]), 
	&(_tmp->get_data_ptr()[batch*prod(matrix_size_os)*domain_size_coils]), &(samples->get_data_ptr()[batch*number_of_samples*domain_size_coils]), 
	raw_pointer_cast(&(*tuples_last)[0]), raw_pointer_cast(&(*bucket_begin)[0]), raw_pointer_cast(&(*bucket_end)[0]),
	double_warp_size_power, half(W), reciprocal(W), huintd_to_reald(matrix_size_os) );
  }
  
  CHECK_FOR_CUDA_ERROR();
  
  bool success = image_wrap( _tmp, image );
  
  delete _tmp;
  
  return success;
}

// TODO: GENERALIZE/TEMPLETIZE the wrapping

// Image wrap kernels

template<class NDTYPE> __inline__ __device__ void 
__wrap_image( uint2 co, uint2 matrix_size_os, uint2 half_wrap, uint2 source_image_size, unsigned int source_image_offset, NDTYPE *source_image, NDTYPE *out )
{
  // Wrap edges
  if( co.x < half_wrap.x )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y),source_image_size)+source_image_offset];
  if( co.x >= (matrix_size_os.x-half_wrap.x) )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y),source_image_size)+source_image_offset];
  if( co.y < half_wrap.y )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x, co.y+half_wrap.y+matrix_size_os.y),source_image_size)+source_image_offset];
  if( co.y >= (matrix_size_os.y-half_wrap.y) )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x, co.y+half_wrap.y-matrix_size_os.y),source_image_size)+source_image_offset];

  // Wrap corners
  if( co.x < half_wrap.x && co.y < half_wrap.y )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y),source_image_size)+source_image_offset];
  if( co.x < half_wrap.x && co.y >= (matrix_size_os.y-half_wrap.y) )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y),source_image_size)+source_image_offset];
  if( co.x >= (matrix_size_os.x-half_wrap.x) && co.y < half_wrap.y )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y),source_image_size)+source_image_offset];
  if( co.x >= (matrix_size_os.x-half_wrap.x) && co.y >= (matrix_size_os.y-half_wrap.y) )
    *out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y),source_image_size)+source_image_offset];
}

template<class NDTYPE> __inline__ __device__ void
__wrap_image( uint3 co, uint3 matrix_size_os, uint3 half_wrap, uint3 source_image_size, unsigned int source_image_offset, NDTYPE *source_image, NDTYPE *out )
{

  // !!! THIS CODE IS ONLY VALID WHEN THE THERE IS NO CONVOLUTION IN THE 'Z' DIMENSION (k-t SENSE) !!!

  // Wrap edges
  if( co.x < half_wrap.x )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y, co.z),source_image_size)+source_image_offset];
  if( co.x >= (matrix_size_os.x-half_wrap.x) )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y, co.z),source_image_size)+source_image_offset];
  if( co.y < half_wrap.y )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x, co.y+half_wrap.y+matrix_size_os.y, co.z),source_image_size)+source_image_offset];
  if( co.y >= (matrix_size_os.y-half_wrap.y) )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x, co.y+half_wrap.y-matrix_size_os.y, co.z),source_image_size)+source_image_offset];

  // Wrap corners
  if( co.x < half_wrap.x && co.y < half_wrap.y )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y, co.z),source_image_size)+source_image_offset];
  if( co.x < half_wrap.x && co.y >= (matrix_size_os.y-half_wrap.y) )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y, co.z),source_image_size)+source_image_offset];
  if( co.x >= (matrix_size_os.x-half_wrap.x) && co.y < half_wrap.y )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y, co.z),source_image_size)+source_image_offset];
  if( co.x >= (matrix_size_os.x-half_wrap.x) && co.y >= (matrix_size_os.y-half_wrap.y) )
    *out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y, co.z),source_image_size)+source_image_offset];

}

template<class NDTYPE> __inline__ __device__ void
__wrap_image( uint4 co, uint4 matrix_size_os, uint4 half_wrap, uint4 source_image_size, unsigned int source_image_offset, NDTYPE *source_image, NDTYPE *out )
{
}

template<class UINTd, class REALd, class REAL, class NDTYPE> __global__ void
image_wrap_kernel( UINTd matrix_size_os, UINTd matrix_size_wrap, unsigned int number_of_images, /*bool accumulate, */NDTYPE *source, NDTYPE *target )
{
  const unsigned int image_idx = blockIdx.y*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x;
  const unsigned int image_number = threadIdx.y;
  const unsigned int idx = image_number*gridDim.y*gridDim.x*blockDim.x+image_idx;

  if( image_idx <prod(matrix_size_os) ){

    const UINTd source_image_size = matrix_size_os+matrix_size_wrap;
    const unsigned int source_image_offset = image_number*prod(source_image_size);
    
    const UINTd half_wrap = matrix_size_wrap>>1;
    const UINTd co = idx_to_co( image_idx, matrix_size_os );
    
    NDTYPE out = source[co_to_idx(co+half_wrap,source_image_size)+source_image_offset];
    
    /*
    if( accumulate )
      out += target_image[idx];
    */

    // Wrap image
    __wrap_image<NDTYPE>( co, matrix_size_os, half_wrap, source_image_size, source_image_offset, source, &out ); 
		
    target[idx] = out;
  }
}

template<class UINTd, class REALd, class REAL, class NDTYPE> bool
NFFT_plan<UINTd, REALd, REAL, NDTYPE>::image_wrap( cuNDArray<NDTYPE> *source, cuNDArray<NDTYPE> *target )
{
  
  // private method - no consistency check. We trust in ourselves.
  
  UINTd source_dims; cuNDA_fromVec( source->get_dimensions(), source_dims );
  UINTd target_dims; cuNDA_fromVec( target->get_dimensions(), target_dims );

  unsigned int number_of_images = (source->get_number_of_dimensions()==d+1) ? source->get_size(source->get_number_of_dimensions()-1) : 1;
  
  // Find dimensions of grid/blocks.

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties( &deviceProp, device );

  unsigned int dimBlockX = ((deviceProp.maxThreadsPerBlock/number_of_images)/deviceProp.warpSize)*deviceProp.warpSize;
  unsigned int last_dim_os = cuNDA_toVec(matrix_size_os)[d-1];

  // TODO: recall why we do this ugly code!

  // Make sure that gridDim.x becomes integer valued
  while(((prod(matrix_size_os)/last_dim_os)%dimBlockX) != 0 ){
    if( dimBlockX<deviceProp.warpSize ){
      printf("\nImplementation error: image_wrap: Cannot reduce the block size below the warp size.\n");
      exit(1);
    }
    else{
      dimBlockX -= deviceProp.warpSize;	
    }
  }

  // TODO: fix me for manyD

  if(dimBlockX == 0){
    printf("\nImplementation error: image_wrap: Too many coils reduces the block size below the warp size.\n");
    exit(1);
  }

  unsigned int last_dim_wrap = cuNDA_toVec(matrix_size_wrap)[d-1];

  if( (sizeof(UINTd)!=sizeof(uint2)) && last_dim_wrap != 0 ){
    printf("\nImplementation error: image_wrap: 3D wrapping not yet implemented.\n");
    exit(1);
  }

  dim3 dimBlock( dimBlockX, number_of_images );
  dim3 dimGrid( (prod(matrix_size_os)/last_dim_os)/dimBlock.x, last_dim_os ); // No 'ceil'!!!

  // Invoke kernel
  image_wrap_kernel<UINTd, REALd, REAL, NDTYPE><<<dimGrid, dimBlock>>>( matrix_size_os, matrix_size_wrap, number_of_images, /*accumulate,*/ source->get_data_ptr(), target->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();

  return true;
}	

//
// Template instantion
//

template class NFFT_plan< uint2, float2, float, cuFloatComplex >;
