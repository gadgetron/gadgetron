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

#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "ndarray_vector_td_utilities.hcu"

#include "vector_td_operators.hcu"
#include "vector_td_utilities.hcu"

#include "check_CUDA.h"

// Includes - CUDA
#include <math_constants.h>
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
NFFT_plan<REAL,D>::NFFT_plan( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, typename uintd<D>::Type fixed_dims, REAL W, int device )
{
  // Minimal initialization
  barebones();

  // Setup plan
  setup( matrix_size, matrix_size_os, fixed_dims, W, device );
  
  if( !initialized )
    cout << endl << "Initialization of the plan failed." << endl;
}

template<class REAL, unsigned int D> 
NFFT_plan<REAL,D>::NFFT_plan( NFFT_plan<REAL,D> *plan )
{
  matrix_size = plan->matrix_size;
  matrix_size_os = plan->matrix_size_os;
  matrix_size_wrap = plan->matrix_size_wrap;
  fixed_dims = plan->fixed_dims;
  non_fixed_dims = plan->non_fixed_dims;
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

template<class REAL, unsigned int D> 
NFFT_plan<REAL,D>::~NFFT_plan()
{
  wipe(NFFT_WIPE_ALL);
}

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, typename uintd<D>::Type fixed_dims, REAL W, int _device )
{	
  // Free memory
  wipe(NFFT_WIPE_ALL);

  //
  // Check if the device is valid
  //

  if( _device<0 )
    cudaGetDevice( &device );
  else
    device = _device;

  cudaDeviceProp deviceProp;  
  cudaGetDeviceProperties( &deviceProp, device );

  unsigned int warp_size = deviceProp.warpSize;
  typename uintd<D>::Type vec_warp_size; to_vector_td<unsigned int,D>(vec_warp_size, warp_size);

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
  this->fixed_dims = fixed_dims;

  for( unsigned int i=0; i<D; i++ ){
    
    if( fixed_dims.vec[i] >= 1 ){
      non_fixed_dims.vec[i] = 0;
    }
    else{
      non_fixed_dims.vec[i] = 1;
    }
  }


  REAL W_half = get_half<REAL>()*W;
  vector_td<REAL,D> W_vec; to_vector_td<REAL,D>(W_vec, W_half);

  to_uintd<REAL,D>( matrix_size_wrap, ceil(W_vec));
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

    if( fixed_dims.vec[dim] > 1 )
      fixed_dims.vec[dim] = 1; // enforce fixed_dims as "boolean style"

    if( fixed_dims.vec[dim] ){
      matrix_size_wrap.vec[dim] = 0;
    }
    
    if( !fixed_dims.vec[dim] &&
	((matrix_size_os.vec[dim]/(matrix_size.vec[dim])) != fracX ||
	 (matrix_size_os.vec[dim]%(matrix_size.vec[dim]) != moduX) )){
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

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory, NFFT_prep_mode mode )
{
  if( trajectory->get_device() != device ){
    cout << endl << "NFFT_plan::preprocess: device mismatch." << endl;
    return false;
  }

  wipe(NFFT_WIPE_PREPROCESSING);

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
  device_vector< vector_td<REAL,D> > trajectory_positions_in( device_pointer_cast< vector_td<REAL,D> >(trajectory->get_data_ptr()), 
							    device_pointer_cast< vector_td<REAL,D> >(trajectory->get_data_ptr()+number_of_samples) );
  switch(D){
    
    // Switch to ensure alignment issues
    
  case 2:
    trajectory_positions = (sizeof(REAL)==sizeof(float)) ? 
      (device_vector< vector_td<REAL,D> >*) (new device_vector<floatd2>(number_of_samples)) : 
      (device_vector< vector_td<REAL,D> >*) (new device_vector<doubled2>(number_of_samples));
    break;
  case 4:
    trajectory_positions = (sizeof(REAL)==sizeof(float)) ? 
      (device_vector< vector_td<REAL,D> >*) (new device_vector<floatd4>(number_of_samples)) : 
      (device_vector< vector_td<REAL,D> >*) (new device_vector<doubled4>(number_of_samples));
    break;
  default:
    trajectory_positions = new device_vector< vector_td<REAL,D> >( number_of_samples );
    break;
  }
  
  CHECK_FOR_CUDA_ERROR();

  vector_td<REAL,D> matrix_size_os_real; to_reald<REAL,unsigned int,D>( matrix_size_os_real, matrix_size_os );
  vector_td<REAL,D> matrix_size_os_plus_wrap_real; to_reald<REAL,unsigned int,D>( matrix_size_os_plus_wrap_real, (matrix_size_os+matrix_size_wrap)>>1 );

  // convert input trajectory in [-1/2;1/2] to [0;matrix_size_os]
  thrust::transform( trajectory_positions_in.begin(), trajectory_positions_in.end(), trajectory_positions->begin(), trajectory_scale<REAL,D>(matrix_size_os_real, matrix_size_os_plus_wrap_real) );

  if( mode != NFFT_PREP_FORWARDS ){
    
    // allocate storage for and compute temporary prefix-sum variable (#cells influenced per sample)
    device_vector<unsigned int> c_p_s(number_of_samples);
    device_vector<unsigned int> c_p_s_ps(number_of_samples);
    
    CHECK_FOR_CUDA_ERROR();
    
    REAL half_W = get_half<REAL>()*W;
    thrust::plus<unsigned int> binary_op;
    thrust::transform(trajectory_positions->begin(), trajectory_positions->end(), c_p_s.begin(), compute_num_cells_per_sample<REAL,D>(half_W));
    inclusive_scan( c_p_s.begin(), c_p_s.end(), c_p_s_ps.begin(), binary_op ); // prefix sum
    
    // Build the vector of (grid_idx, sample_idx) tuples. Actually kept in two seperate vectors.
    unsigned int num_pairs = c_p_s_ps.back();
    thrust::device_vector<unsigned int> tuples_first = device_vector<unsigned int>(num_pairs);
    tuples_last = new device_vector<unsigned int>(num_pairs);
    
    CHECK_FOR_CUDA_ERROR();
    
    // Fill tuple vector
    write_pairs<REAL,D>( matrix_size_os, matrix_size_wrap, number_of_samples, W, 
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

  if( mode != NFFT_PREP_FORWARDS )
    preprocessed_NFFT_H = true;

  return true;
}


template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::compute( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode )
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
  cuNDArray<typename complext<REAL>::Type> *working_image = 0x0, *working_samples = 0x0;
  
  typename uintd<D>::Type image_dims; 
  if( !cuNDA_fromVec<D>( image->get_dimensions(), image_dims ) ){
    cout << "NFFT_plan::compute: image dimension undermatch the plan" << endl;
    return false;
  }
  
  bool oversampled_image = (image_dims==matrix_size_os); 
  unsigned int num_batches = (image->get_number_of_dimensions()>D) ? image->get_size(image->get_number_of_dimensions()-1) : 1;
  vector<unsigned int> vec_dims = cuNDA_toVec<D>(matrix_size_os); 
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

      working_image = cuNDArray<typename complext<REAL>::Type>::allocate(vec_dims);

      cuNDA_expand_with_zero_fill<typename complext<REAL>::Type,D>( image, working_image );
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
      working_samples = new cuNDArray<typename complext<REAL>::Type>(*samples); 
   
      cuNDA_scale( weights, working_samples );
    }
    else{
      working_samples = samples;
    }
    
    if( !oversampled_image ){

      working_image = cuNDArray<typename complext<REAL>::Type>::allocate(vec_dims); 
    }
    else{
      working_image = image;
    }

    success = compute_NFFT_H( working_samples, working_image );

    if( success ){
      if( !oversampled_image ){
	cuNDA_crop<typename complext<REAL>::Type,D>( (matrix_size_os-matrix_size)>>1, working_image, image );
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

  CHECK_FOR_CUDA_ERROR();
  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::compute_iteration( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode )
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
  cuNDArray<typename complext<REAL>::Type> *working_image = 0x0, *working_samples = 0x0;

  typename uintd<D>::Type image_dims; 
  if( !cuNDA_fromVec<D>( image->get_dimensions(), image_dims ) ){
    cout << "NFFT_plan::compute_iteration: image dimension undermatch the plan" << endl;
    return false;
  }
  
  bool oversampled_image = (image_dims==matrix_size_os); 
  unsigned int num_batches = (image->get_number_of_dimensions()>D) ? image->get_size(image->get_number_of_dimensions()-1) : 1;
  vector<unsigned int> vec_dims = cuNDA_toVec<D>(matrix_size_os); 
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
      
      working_image = cuNDArray<typename complext<REAL>::Type>::allocate(vec_dims); 

      cuNDA_expand_with_zero_fill<typename complext<REAL>::Type,D>( image, working_image );
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
	cuNDA_crop<typename complext<REAL>::Type,D>( (matrix_size_os-matrix_size)>>1, working_image, image );
	delete working_image; working_image = 0x0;
      }
    }
    
    break;

  case NFFT_BACKWARDS: // iteration from samples

    // Density compensation
    if( weights ){
      working_samples = new cuNDArray<typename complext<REAL>::Type>(*samples);
      cuNDA_scale( weights, working_samples );
    }
    else{
      working_samples = samples;
    }
    
    if( !oversampled_image )
      
      working_image = cuNDArray<typename complext<REAL>::Type>::allocate(vec_dims); 

    else
      working_image = image;

    success = compute_NFFT_H( working_samples, working_image );
    
    if( success ){
      
      cuNDA_zero_fill_border<typename complext<REAL>::Type,D>( matrix_size, working_image );
      
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

  CHECK_FOR_CUDA_ERROR();
  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::convolve( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image, cuNDArray<REAL> *weights, NFFT_comp_mode mode )
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
  cuNDArray<typename complext<REAL>::Type> *working_image = 0x0, *working_samples = 0x0;
  
  typename uintd<D>::Type image_dims; 
  if( !cuNDA_fromVec<D>( image->get_dimensions(), image_dims ) ){
    cout << "NFFT_plan::convolve: image dimension undermatch the plan" << endl;
    return false;
  }
  
  bool oversampled_image = (image_dims==matrix_size_os); 
  unsigned int num_batches = (image->get_number_of_dimensions()>D) ? image->get_size(image->get_number_of_dimensions()-1) : 1;
  vector<unsigned int> vec_dims = cuNDA_toVec<D>(matrix_size_os); 
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

      working_image = cuNDArray<typename complext<REAL>::Type>::allocate(vec_dims); 

      cuNDA_expand_with_zero_fill<typename complext<REAL>::Type,D>( image, working_image );
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
      working_samples = new cuNDArray<typename complext<REAL>::Type>(*samples);
      cuNDA_scale( weights, working_samples );
    }
    else{
      working_samples = samples;
    }
    
    if( !oversampled_image ){

      working_image = cuNDArray<typename complext<REAL>::Type>::allocate(vec_dims); 
    }
    else{
      working_image = image;
    }
    
    success = convolve_NFFT_H( working_samples, working_image );

    if( success ){
      if( !oversampled_image ){
	cuNDA_crop<typename complext<REAL>::Type,D>( (matrix_size_os-matrix_size)>>1, working_image, image );
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

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::FFT(cuNDArray<typename complext<REAL>::Type> *data, NFFT_comp_mode mode, bool do_scale )
{
  if( data->get_device() != device ){
    cout << endl << "NFFT_plan::FFT: device mismatch." << endl;
    return false;
  }

  typename uintd<D>::Type _dims_to_transform = counting_vec<D>();
  vector<unsigned int> dims_to_transform = cuNDA_toVec<D>( _dims_to_transform );

  int res;
  if( mode == NFFT_FORWARDS ){
    res = cuNDFFT().fft( (cuNDArray<cuFloatComplex>*)data, dims_to_transform ); // TODO: remove casting / fix fft interface in cuNDFFT.h/cu
  }
  else{
    res = cuNDFFT().ifft( (cuNDArray<cuFloatComplex>*)data, dims_to_transform, do_scale );
  }

  if( res == 0 )
    return true;
  else
    return false;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::deapodize( cuNDArray<typename complext<REAL>::Type> *image )
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

  if( !((image->get_number_of_dimensions() == D) || (image->get_number_of_dimensions() == (D+1))) ){
    cout << endl << "NFFT_plan: Number of image dimensions mismatch the plan." << endl;
    return false;
  }    

  typename uintd<D>::Type image_dims; cuNDA_fromVec<D>( image->get_dimensions(), image_dims );
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
    
    if( samples->get_number_of_dimensions()==2 || image->get_number_of_dimensions() == D+1 ){
      if( !(samples->get_number_of_dimensions()==2) || !(image->get_number_of_dimensions() == D+1) || !(samples->get_size(1)==image->get_size(D)) ){
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

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::barebones()
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

template<class REAL, unsigned int D> 
void NFFT_plan<REAL,D>::wipe( NFFT_wipe_mode mode )
{
  if( mode==NFFT_WIPE_ALL && initialized ){
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

template<class REAL, unsigned int D> 
bool NFFT_plan<REAL,D>::compute_beta()
{	
  // Compute Kaiser-Bessel beta paramter according to the formula provided in 
  // Beatty et. al. IEEE TMI 2005;24(6):799-808.

  beta = (M_PI*sqrt((W*W)/(alpha*alpha)*(alpha-0.5)*(alpha-0.5)-0.8)); 

  return true;
}

//
// Grid fictitious trajectory with a single sample at the origin
//

template<class REAL, unsigned int D> __global__ void
compute_deapodization_filter_kernel( typename uintd<D>::Type matrix_size_os, typename reald<REAL,D>::Type matrix_size_os_real, 
				     typename uintd<D>::Type fixed_dims, 
				     REAL W, REAL half_W, REAL one_over_W, REAL beta, typename complext<REAL>::Type*image_os )
{

  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int num_elements = prod(matrix_size_os);

  if( idx <num_elements ){

    // Compute weight from Kaiser-Bessel filter
    const vector_td<unsigned int,D> cell_pos = idx_to_co<D>(idx, matrix_size_os);

    // Sample position ("origin")
    const vector_td<REAL,D> sample_pos = get_half<REAL>()*matrix_size_os_real;

    // Calculate the distance between the cell and the sample
    vector_td<REAL,D> cell_pos_real; to_reald<REAL,unsigned int,D>(cell_pos_real, cell_pos);
    const typename reald<REAL,D>::Type delta = abs(sample_pos-cell_pos_real);

    // Compute convolution weight. 
    REAL weight; 
    REAL zero = get_zero<REAL>();
    vector_td<REAL,D> half_W_vec; to_vector_td<REAL,D>( half_W_vec, half_W );

    if( weak_greater(delta, half_W_vec ) )
      weight = zero;
    else{ 
      weight = KaiserBessel<REAL>( delta, matrix_size_os_real, one_over_W, beta, fixed_dims );
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

  if( initialized && deapodization_filter ) 
    delete( deapodization_filter );
  
  deapodization_filter = cuNDArray<typename complext<REAL>::Type>::allocate(cuNDA_toVec<D>(matrix_size_os)); 
  
  if( !deapodization_filter ){
    cout << endl << "cuNDArray allocation failed" << endl;
    return false;
  }

  vector_td<REAL,D> matrix_size_os_real; to_reald<REAL,unsigned int, D>(matrix_size_os_real, matrix_size_os);
  
  // Be optimistic
  bool success = true;

  // Find dimensions of grid/blocks.

  dim3 dimBlock( 256 );
  dim3 dimGrid( (unsigned int) ceil((double)prod(matrix_size_os)/(double)dimBlock.x) );

  // Invoke kernel
  compute_deapodization_filter_kernel<REAL,D><<<dimGrid, dimBlock>>> 
    ( matrix_size_os, matrix_size_os_real, fixed_dims, W, get_half<REAL>()*W, reciprocal<REAL>(W), beta, deapodization_filter->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();
  
  // FFT
  if( success )
    success = FFT( deapodization_filter, NFFT_BACKWARDS, false );
  
  // Scale (multiplication by N, i.e. not what the FFT provides)
  // Now the NFFT achieves scaling by applying the deapodization filter and we save a little bit of time...
  if( success )
    cuNDA_scale<REAL,typename complext<REAL>::Type>( (REAL)prod(matrix_size_os), deapodization_filter );
  
  // Reciprocal
  if( success )
    cuNDA_reciprocal<typename complext<REAL>::Type>( deapodization_filter );

  if( !success ){
    delete deapodization_filter;
  }
  
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
    success = FFT( image, NFFT_FORWARDS ); 

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
    success = FFT( image, NFFT_BACKWARDS, false ); // scaling is cared for by deapodization
  
  // Deapodization
  if( success )
    success = deapodize( image );
  
  return success;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::convolve_NFFT( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image )
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

  switch(D){
    
  case 2:
    if( generation == 1 ){
      threads_per_block = NFFT_THREADS_PER_2D_KERNEL_1x;
      max_coils = NFFT_MAX_COILS_2D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_THREADS_PER_2D_KERNEL_2x;
      max_coils = NFFT_MAX_COILS_2D_COMPUTE_2x;
    }
    break;
    
  case 3:
    if( generation == 1 ){
      threads_per_block = NFFT_THREADS_PER_2D_KERNEL_1x;
      max_coils = NFFT_MAX_COILS_2D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_THREADS_PER_2D_KERNEL_2x;
      max_coils = NFFT_MAX_COILS_2D_COMPUTE_2x;
    }
    break;
    
  case 4:
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
  
  // We can (only) convolve max_coils batches per run due to shared memory issues. 
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
  size_t bytes_per_thread = domain_size_coils * sizeof( vector_td<REAL,D> );
  size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof( vector_td<REAL,D> );

  /*
    Invoke kernel
  */

  unsigned int double_warp_size_power=0;
  unsigned int __tmp = warp_size<<1;
  while(__tmp!=1){
    __tmp>>=1;
    double_warp_size_power++;
  }
  
  vector_td<REAL,D> matrix_size_os_real; to_reald<REAL,unsigned int,D>( matrix_size_os_real, matrix_size_os );

  for( unsigned int batch = 0; batch<num_repetitions; batch++ ){
    NFFT_convolve_kernel<REAL,D>
      <<<dimGrid, dimBlock, (batch==num_repetitions-1) ? dimBlock.x*bytes_per_thread_tail : dimBlock.x*bytes_per_thread>>>
      ( alpha, beta, W, matrix_size_os, matrix_size_wrap, fixed_dims, number_of_samples, 
	(batch==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
	raw_pointer_cast(&(*trajectory_positions)[0]), 
	&(image->get_data_ptr()[batch*prod(matrix_size_os)*domain_size_coils]), &(samples->get_data_ptr()[batch*number_of_samples*domain_size_coils]), 
	double_warp_size_power, get_half<REAL>()*W, reciprocal(W), matrix_size_os_real, non_fixed_dims );
  }
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::convolve_NFFT_H( cuNDArray<typename complext<REAL>::Type> *samples, cuNDArray<typename complext<REAL>::Type> *image )
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

  switch(D){

  case 2:
    if( generation == 1 ){
      threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_1x;
      max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_2x;
      max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_2x;
    }
    break;
      
  case 3:
    if( generation == 1 ){
      threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_1x;
      max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_1x;
    }
    else{
      threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_2x;
      max_coils = NFFT_H_MAX_COILS_2D_COMPUTE_2x;
    }
    break;
    
  case 4:
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
  size_t bytes_per_thread = domain_size_coils * sizeof( vector_td<REAL,D> );
  size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof( vector_td<REAL,D> );

  vector<unsigned int> vec_dims = cuNDA_toVec<D>(matrix_size_os+matrix_size_wrap); 
  if( domain_size_coils_desired > 1 ) 
    vec_dims.push_back(domain_size_coils_desired);
  
  // Allocate and clear temporary image that includes a wrapping zone
  cuNDArray<typename complext<REAL>::Type> *_tmp = cuNDArray<typename complext<REAL>::Type>::allocate(vec_dims); 
  
  if( !_tmp ){
    cout << endl << "Memory allocation failed before convolution." << endl;
    return false;
  }
  
  unsigned int double_warp_size_power=0, __tmp = warp_size<<1;
  while(__tmp!=1){
    __tmp>>=1;
    double_warp_size_power++;
  }

  vector_td<REAL,D> matrix_size_os_real; to_reald<REAL,unsigned int,D>( matrix_size_os_real, matrix_size_os );

  for( unsigned int batch = 0; batch<num_repetitions; batch++ ){
    NFFT_H_convolve_kernel<REAL,D>
      <<<dimGrid, dimBlock, (batch==num_repetitions-1) ? dimBlock.x*bytes_per_thread_tail : dimBlock.x*bytes_per_thread>>>
      ( alpha, beta, W, matrix_size_os+matrix_size_wrap, fixed_dims, number_of_samples, 
	(batch==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
	raw_pointer_cast(&(*trajectory_positions)[0]), 
	&(_tmp->get_data_ptr()[batch*prod(matrix_size_os)*domain_size_coils]), &(samples->get_data_ptr()[batch*number_of_samples*domain_size_coils]), 
	raw_pointer_cast(&(*tuples_last)[0]), raw_pointer_cast(&(*bucket_begin)[0]), raw_pointer_cast(&(*bucket_end)[0]),
	double_warp_size_power, get_half<REAL>()*W, reciprocal(W), matrix_size_os_real );
  }
  
  //cudaThreadSynchronize();
  CHECK_FOR_CUDA_ERROR();
  
  bool success = image_wrap( _tmp, image );
 
  delete _tmp;

  return success;
}

// Image wrap kernels

template<class REAL, unsigned int D> __global__ void
image_wrap_kernel( vector_td<unsigned int,D> matrix_size_os, vector_td<unsigned int,D> non_fixed_dims, vector_td<unsigned int,D> matrix_size_wrap, 
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

  const unsigned int image_offset_src = threadIdx.z * num_elements_per_image_src;
  
  const vector_td<unsigned int,D> co = idx_to_co<D>(idx, matrix_size_os);
  const vector_td<unsigned int,D> half_wrap = matrix_size_wrap>>1;
  
  // Make "boolean" vectors denoting whether wrapping needs to be performed in a given direction (forwards/backwards)
  vector_td<unsigned int,D> B_l = vector_less( co, half_wrap ); 
  vector_td<unsigned int,D> B_r = vector_greater_equal( co, matrix_size_os-half_wrap ); 
  B_l *= non_fixed_dims; B_r *= non_fixed_dims; // don't wrap fixed dimensions
  
  typename complext<REAL>::Type result = in[co_to_idx<D>(co+half_wrap, matrix_size_os+matrix_size_wrap) + image_offset_src];

  if( sum(B_l+B_r) > 0 ){
    
    // Fold back the wrapping zone onto the image ("periodically")
    //
    // There is 2^D-1 ways to pick combinations of dimensions in D-dimensionsal space, e.g. 
    // 
    //  { x, y, xy } in 2D
    //  { x, y, x, xy, xz, yz, xyz } in 3D
    //
    // Every "letter" in aech combination provides two possible wraps (eiher end of the dimension)
    // 
    // For every 2^D-1 combinations DO
    //   - find the number of dimensions, d, in the combination
    //   - create 2^(d) stride vectors and test for wrapping the 'B'-vectors above.
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
	  vector_td<int,D> src_co_int; to_intd(src_co_int, co+half_wrap);
	  vector_td<int,D> matrix_size_os_int; to_intd(matrix_size_os_int, matrix_size_os);
	  vector_td<int,D> co_offset_int = src_co_int + stride*matrix_size_os_int;
	  vector_td<unsigned int,D> co_offset; to_uintd(co_offset, co_offset_int);
	  result += in[co_to_idx<D>(co_offset, matrix_size_os+matrix_size_wrap) + image_offset_src];
	  break; // only one stride per combination can contribute (e.g. one edge, one corner)
	} 
      } 
    }
  }
  
  // Output
  const unsigned int image_offset_tgt = threadIdx.z * prod(matrix_size_os);
  out[idx+image_offset_tgt] = result;
}

template<class REAL, unsigned int D> bool
NFFT_plan<REAL,D>::image_wrap( cuNDArray<typename complext<REAL>::Type> *source, cuNDArray<typename complext<REAL>::Type> *target )
{
  unsigned int number_of_images = (source->get_number_of_dimensions()==(D+1)) ? source->get_size(source->get_number_of_dimensions()-1) : 1;
  
  // Set dimensions of grid/blocks.
  unsigned int bdim = 16;
  dim3 dimBlock( bdim, (number_of_images==1) ? bdim : 1, number_of_images );
  dim3 dimGrid( prod(matrix_size_os)/(dimBlock.x*dimBlock.y) );

  // Safety check
  if( (matrix_size_os.vec[0]%bdim) || (number_of_images==1 && (matrix_size_os.vec[1]%bdim)) ) {
    cout << endl << "dimensions mismatch in wrapping setup." << endl;
    return false;
  }

  // Invoke kernel
  image_wrap_kernel<REAL,D><<<dimGrid, dimBlock>>>( matrix_size_os, non_fixed_dims, matrix_size_wrap, /*accumulate,*/
						    source->get_data_ptr(), target->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();

  return true;
}	

//
// Template instantion
//

template class NFFT_plan< float, 2 >;
template class NFFT_plan< double, 2 >;

template class NFFT_plan< float, 3 >;
template class NFFT_plan< double, 3 >;

template class NFFT_plan< float, 4 >;
template class NFFT_plan< double, 4 >;
