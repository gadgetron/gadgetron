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

// Includes - Thrust
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/extrema.h>
// Includes - Gadgetron
#include "cuNFFT.h"
#include "cuNDFFT.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "vector_td_io.h"
#include "cudaDeviceManager.h"
#include "check_CUDA.h"

// Includes - CUDA
#include <device_functions.h>
#include <math_constants.h>
#include <cufft.h>


// Includes - stdlibs
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <cmath>
#include <sstream>
#include <stdexcept>

//using namespace std;
using std::vector;
using namespace thrust;
using namespace Gadgetron;

// Kernel configuration  
#define NFFT_MAX_COILS_COMPUTE_1x    8
#define NFFT_MAX_COILS_COMPUTE_2x   16
#define NFFT_THREADS_PER_KERNEL    192

// Reference to shared memory
extern __shared__ char _shared_mem[];

// Includes containing the NFFT convolution implementation
#include "KaiserBessel_kernel.cu"
#include "NFFT_C2NC_conv_kernel.cu"
#include "NFFT_NC2C_conv_kernel.cu"
#include "NFFT_NC2C_atomic_conv_kernel.cu"
#include "NFFT_preprocess_kernel.cu"

// Default template arguments requires c++-0x ?
typedef float dummy;

// The declaration of atomic/non-atomic NC2C convolution
// We would love to hide this inside the class, but the compiler core dumps on us when we try...
//
template<class REAL, unsigned int D, bool ATOMICS> struct _convolve_NFFT_NC2C{
  static bool apply( cuNFFT_plan<REAL,D,ATOMICS> *plan, 
                     cuNDArray<complext<REAL> > *in, 
                     cuNDArray<complext<REAL> > *out, 
                     bool accumulate );
};

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
    throw cuda_error("Error: cuNFFT : unable to get device no");
  }

  if( device != *old_device && cudaSetDevice(device) != cudaSuccess) {
    throw cuda_error("Error : cuNFFT : unable to set device no");
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
    *out = *in1_int; // device transfer by assignment
  } 
  
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
    throw cuda_error("Error: cuNFFT : unable to get device no");
  }
  
  // Restore old device
  if( device != old_device && cudaSetDevice(old_device) != cudaSuccess) {
    throw cuda_error("Error: cuNFFT : unable to restore device no");
  }
  
  return true;
}


//
// Public class methods
//

template<class REAL, unsigned int D, bool ATOMICS> 
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::cuNFFT_plan()
{
  // Minimal initialization
  barebones();
}

template<class REAL, unsigned int D, bool ATOMICS> 
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::cuNFFT_plan( typename uint64d<D>::Type matrix_size, typename uint64d<D>::Type matrix_size_os, REAL W, int device )
{
  // Minimal initialization
  barebones();

  // Setup plan
  setup( matrix_size, matrix_size_os, W, device );
}

template<class REAL, unsigned int D, bool ATOMICS> 
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::~cuNFFT_plan()
{
  wipe(NFFT_WIPE_ALL);
}

template<class REAL, unsigned int D, bool ATOMICS> 
void Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::setup( typename uint64d<D>::Type matrix_size, typename uint64d<D>::Type matrix_size_os, REAL W, int _device )
{
  // Free memory
  wipe(NFFT_WIPE_ALL);

  //
  // Check if the device is valid
  //

  if( _device<0 ){
    if( cudaGetDevice( &device ) != cudaSuccess ){
      throw cuda_error("Error: cuNFFT_plan::setup: unable to determine device properties.");
    }
  }
  else
    device = _device;

  // The convolution does not work properly for very small convolution kernel widths
  // (experimentally observed limit)

  if( W < REAL(1.8) ) {
    throw std::runtime_error("Error: the convolution kernel width for the cuNFFT plan is too small.");
  }

  typename uint64d<D>::Type vec_warp_size( (size_t)(cudaDeviceManager::Instance()->warp_size(device)) );

  //
  // Check input against certain requirements
  //
  
  if( sum(matrix_size%vec_warp_size) || sum(matrix_size_os%vec_warp_size) ){
    //GDEBUG_STREAM("Matrix size: " << matrix_size << std::endl);
    //GDEBUG_STREAM("Matrix size os: " << matrix_size_os << std::endl);
    //GDEBUG_STREAM("Warp size: " << vec_warp_size << std::endl);
    throw std::runtime_error("Error: Illegal matrix size for the cuNFFT plan (not a multiple of the warp size)");
  }

  //
  // Setup private variables
  //

  this->matrix_size = matrix_size;
  this->matrix_size_os = matrix_size_os;

  REAL W_half = REAL(0.5)*W;
  vector_td<REAL,D> W_vec(W_half);

  matrix_size_wrap = vector_td<size_t,D>( ceil(W_vec) );
  matrix_size_wrap<<=1; 
  
  alpha = vector_td<REAL,D>(matrix_size_os) / vector_td<REAL,D>(matrix_size);
  
  typename reald<REAL,D>::Type ones(REAL(1));
  if( weak_less( alpha, ones ) ){
    throw std::runtime_error("Error: cuNFFT : Illegal oversampling ratio suggested");
  }

  this->W = W;
  
  // Compute Kaiser-Bessel beta
  compute_beta();
  
  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    throw cuda_error("Error: cuNFFT_plan::setup: unable to get device no");
  }  
  if( device != device_no_old && cudaSetDevice(device) != cudaSuccess) {
    throw cuda_error("Error: cuNFFT_plan::setup: unable to set device");
  }  
  initialized = true;

  if( device != device_no_old && cudaSetDevice(device_no_old) != cudaSuccess) {
    throw cuda_error("Error: cuNFFT_plan::setup: unable to restore device");
  }
}

template<class REAL, unsigned int D, bool ATOMICS> 
void Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory, NFFT_prep_mode mode )
{
  if( !trajectory || trajectory->get_number_of_elements()==0 ){
    throw std::runtime_error("Error: cuNFFT_plan::preprocess: invalid trajectory");
  }
  
  if( !initialized ){
    throw std::runtime_error("Error: cuNFFT_plan::preprocess: cuNFFT_plan::setup must be invoked prior to preprocessing.");
  }

  wipe(NFFT_WIPE_PREPROCESSING);

  cuNDArray<typename reald<REAL,D>::Type> *trajectory_int;
  int old_device;

  if( !prepare<typename reald<REAL,D>::Type,dummy,dummy>(device, &old_device, trajectory, &trajectory_int ) ){
    throw cuda_error("Error: cuNFFT_plan::preprocess: device preparation error.");
  }
    
  number_of_samples = trajectory_int->get_size(0);
  number_of_frames = trajectory_int->get_number_of_elements()/number_of_samples;

  // Make sure that the trajectory values are within range [-1/2;1/2]
  thrust::pair< thrust::device_ptr<REAL>, thrust::device_ptr<REAL> > mm_pair = 
    thrust::minmax_element( device_pointer_cast<REAL>((REAL*)trajectory_int->get_data_ptr()), 
                            device_pointer_cast<REAL>(((REAL*)trajectory_int->get_data_ptr())+trajectory_int->get_number_of_elements()*D ));
  
  if( *mm_pair.first < REAL(-0.5) || *mm_pair.second > REAL(0.5) ){
	  std::stringstream ss;
	  ss << "Error: cuNFFT::preprocess : trajectory [" << *mm_pair.first << "; " << *mm_pair.second << "] out of range [-1/2;1/2]";
    throw std::runtime_error(ss.str());
  }
  
  // Make Thrust device vector of trajectory and samples
  device_vector< vector_td<REAL,D> > trajectory_positions_in
    ( device_pointer_cast< vector_td<REAL,D> >(trajectory_int->get_data_ptr()), 
      device_pointer_cast< vector_td<REAL,D> >(trajectory_int->get_data_ptr()+trajectory_int->get_number_of_elements() ));
  
  trajectory_positions = new device_vector< vector_td<REAL,D> >( trajectory_int->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  vector_td<REAL,D> matrix_size_os_real = vector_td<REAL,D>( matrix_size_os );
  vector_td<REAL,D> matrix_size_os_plus_wrap_real = vector_td<REAL,D>( (matrix_size_os+matrix_size_wrap)>>1 );

  // convert input trajectory in [-1/2;1/2] to [0;matrix_size_os]
  thrust::transform( trajectory_positions_in.begin(), trajectory_positions_in.end(), trajectory_positions->begin(), 
                     trajectory_scale<REAL,D>(matrix_size_os_real, matrix_size_os_plus_wrap_real) );
  
  CHECK_FOR_CUDA_ERROR();

  if( !( mode == NFFT_PREP_C2NC || ATOMICS )){

    // allocate storage for and compute temporary prefix-sum variable (#cells influenced per sample)
    device_vector<unsigned int> c_p_s(trajectory_int->get_number_of_elements());
    device_vector<unsigned int> c_p_s_ps(trajectory_int->get_number_of_elements());
    CHECK_FOR_CUDA_ERROR();
    
    REAL half_W = REAL(0.5)*W;
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
    write_pairs<REAL,D>( vector_td<unsigned int,D>(matrix_size_os), vector_td<unsigned int,D>(matrix_size_wrap), number_of_samples, number_of_frames, W,
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

  preprocessed_C2NC = true;

  if( mode != NFFT_PREP_C2NC )
    preprocessed_NC2C = true;

  if( !restore<typename reald<REAL,D>::Type,dummy,dummy>(old_device, trajectory, trajectory, trajectory_int) ){
    throw cuda_error("Error: cuNFFT_plan::preprocess: unable to restore compute device.");
  }
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::compute( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out,
                                                 cuNDArray<REAL> *dcw, NFFT_comp_mode mode )
{  
  // Validity checks
  
  unsigned char components;

  if( mode == NFFT_FORWARDS_C2NC ) 
    components = _NFFT_CONV_C2NC + _NFFT_FFT + _NFFT_DEAPODIZATION;

  else if( mode == NFFT_FORWARDS_NC2C ) 
    components = _NFFT_CONV_NC2C + _NFFT_FFT + _NFFT_DEAPODIZATION;

  else if( mode == NFFT_BACKWARDS_NC2C ) 
    components = _NFFT_CONV_NC2C + _NFFT_FFT + _NFFT_DEAPODIZATION;

  else if( mode == NFFT_BACKWARDS_C2NC ) 
    components = _NFFT_CONV_C2NC + _NFFT_FFT + _NFFT_DEAPODIZATION;
  else{
    throw std::runtime_error("Error: cuNFFT_plan::compute: unknown mode");
  }
  
  {
    cuNDArray<complext<REAL> > *samples, *image;

    if( mode == NFFT_FORWARDS_C2NC || mode == NFFT_BACKWARDS_C2NC ){
      image = in; samples = out;
    } else{
      image = out; samples = in;
    }
    
    check_consistency( samples, image, dcw, components );
  }
  
  cuNDArray<complext<REAL> > *in_int = 0x0, *out_int = 0x0;
  cuNDArray<REAL> *dcw_int = 0x0;
  int old_device;

  if( !prepare<complext<REAL>, complext<REAL>, REAL>
      (device, &old_device, in, &in_int, out, &out_int, dcw, &dcw_int ) ){
    throw cuda_error("Error: cuNFFT_plan::compute: device preparation error.");
  }

  typename uint64d<D>::Type image_dims = from_std_vector<size_t,D>
    ( (mode == NFFT_FORWARDS_C2NC || mode == NFFT_BACKWARDS_C2NC ) ? *in->get_dimensions() : *out->get_dimensions() );
  bool oversampled_image = (image_dims==matrix_size_os);
  
  vector<size_t> vec_dims = to_std_vector(matrix_size_os);
  {
    cuNDArray<complext<REAL> > *image = ((mode == NFFT_FORWARDS_C2NC || mode == NFFT_BACKWARDS_C2NC ) ? in : out );
    for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
      vec_dims.push_back(image->get_size(d));
  }

  cuNDArray<complext<REAL> > *working_image = 0x0, *working_samples = 0x0;

  switch(mode){

  case NFFT_FORWARDS_C2NC:
    
    if( !oversampled_image ){
      working_image = new cuNDArray<complext<REAL> >(&vec_dims);
      pad<complext<REAL>, D>( in_int, working_image );
    }
    else{
      working_image = in_int;
    }
    
    compute_NFFT_C2NC( working_image, out_int );

    if( dcw_int )
        	*out_int *= *dcw_int;

    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }    
    break;
    
  case NFFT_FORWARDS_NC2C:

    // Density compensation
    if( dcw_int ){
      working_samples = new cuNDArray<complext<REAL> >(*in_int);
      *working_samples *= *dcw_int;
    }
    else{
      working_samples = in_int;
    }
    
    if( !oversampled_image ){
      working_image = new cuNDArray<complext<REAL> >(&vec_dims);
    }
    else{
      working_image = out_int;
    }

    compute_NFFT_NC2C( working_samples, working_image );

    if( !oversampled_image ){
      crop<complext<REAL>, D>( (matrix_size_os-matrix_size)>>1, working_image, out_int );
    }
    
    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    if( dcw_int ){
      delete working_samples; working_samples = 0x0;
    }    
    break;
    
  case NFFT_BACKWARDS_NC2C:
    
    // Density compensation
    if( dcw_int ){
      working_samples = new cuNDArray<complext<REAL> >(*in_int);
      *working_samples *= *dcw_int;
    }
    else{
      working_samples = in_int;
    }
    
    if( !oversampled_image ){
      working_image = new cuNDArray<complext<REAL> >(&vec_dims);
    }
    else{
      working_image = out_int;
    }
    
    compute_NFFTH_NC2C( working_samples, working_image );
    
    if( !oversampled_image ){
      crop<complext<REAL> ,D>( (matrix_size_os-matrix_size)>>1, working_image, out_int );
    }
    
    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    if( dcw_int ){
      delete working_samples; working_samples = 0x0;
    }    
    break;
    
  case NFFT_BACKWARDS_C2NC:
    
    if( !oversampled_image ){
      working_image = new cuNDArray<complext<REAL> >(&vec_dims);
      
      pad<complext<REAL>, D>( in_int, working_image );
    }
    else{
      working_image = in_int;
    }
    
    compute_NFFTH_C2NC( working_image, out_int );
    
    if( dcw_int )
    	*out_int *= *dcw_int;



    if( !oversampled_image ){
      delete working_image; working_image = 0x0;
    }
    
    break;
  };
  
  if( !restore<complext<REAL> ,complext<REAL> ,REAL>
      (old_device, out, out, out_int, in, in_int, dcw, dcw_int ) ){
    throw cuda_error("Error: cuNFFT_plan::compute: unable to restore compute device.");
  }
  
  CHECK_FOR_CUDA_ERROR();
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::mult_MH_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out,
                                                   cuNDArray<REAL> *dcw, std::vector<size_t> halfway_dims )
{
  // Validity checks
  
  unsigned char components = _NFFT_CONV_C2NC + _NFFT_CONV_NC2C + _NFFT_FFT + _NFFT_DEAPODIZATION;
  
  if( in->get_number_of_elements() != out->get_number_of_elements() ){
    throw std::runtime_error("Error: cuNFFT_plan::mult_MH_M: in/out image sizes mismatch");
  }
  
  cuNDArray<complext<REAL> > *working_samples = new cuNDArray<complext<REAL> >(&halfway_dims);

  check_consistency( working_samples, in, dcw, components );
  
  cuNDArray<complext<REAL> > *in_int = 0x0;
  cuNDArray<complext<REAL> > *out_int = 0x0;
  cuNDArray<REAL> *dcw_int = 0x0;
  int old_device;
  
  if( !prepare<complext<REAL>, complext<REAL>, REAL>
      (device, &old_device, in, &in_int, out, &out_int, dcw, &dcw_int ) ){
    throw cuda_error("Error: cuNFFT_plan::mult_MH_M: device preparation error.");
  }
  
  cuNDArray<complext<REAL> > *working_image = 0x0;

  typename uint64d<D>::Type image_dims = from_std_vector<size_t,D>(*in->get_dimensions()); 
  bool oversampled_image = (image_dims==matrix_size_os); 
 
  vector<size_t> vec_dims = to_std_vector(matrix_size_os); 
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ )
    vec_dims.push_back(in->get_size(d));
  
  if( !oversampled_image ){
    working_image = new cuNDArray<complext<REAL> >(&vec_dims);
    pad<complext<REAL>, D>( in_int, working_image );
  }
  else{
    working_image = in_int;
  }
  
  compute_NFFT_C2NC( working_image, working_samples );
  
  // Density compensation
  if( dcw ){
    *working_samples *= *dcw_int;
    *working_samples *= *dcw_int;
  }
    
  compute_NFFTH_NC2C( working_samples, working_image );
    
  delete working_samples;
  working_samples = 0x0;
    
  if( !oversampled_image ){
    crop<complext<REAL>, D>( (matrix_size_os-matrix_size)>>1, working_image, out_int );
    delete working_image; working_image = 0x0;
  }
        
  restore<complext<REAL> ,complext<REAL> ,REAL>
    (old_device, out, out, out_int, in, in_int, dcw, dcw_int );
    
  CHECK_FOR_CUDA_ERROR();
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::convolve( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out,
                                                  cuNDArray<REAL> *dcw, NFFT_conv_mode mode, bool accumulate )
{
  unsigned char components;

  if( mode == NFFT_CONV_C2NC ) 
    components = _NFFT_CONV_C2NC;
  else
    components = _NFFT_CONV_NC2C;
  
  {
    cuNDArray<complext<REAL> > *samples, *image;
    
    if( mode == NFFT_CONV_C2NC ){
      image = in; samples = out;
    } else{
      image = out; samples = in;
    }
    
    check_consistency( samples, image, dcw, components );
  }
  
  cuNDArray<complext<REAL> > *in_int = 0x0, *out_int = 0x0;
  cuNDArray<REAL> *dcw_int = 0x0;
  int old_device;
  
  prepare<complext<REAL>, complext<REAL>, REAL>
    (device, &old_device, in, &in_int, out, &out_int, dcw, &dcw_int );
  
  cuNDArray<complext<REAL> > *working_samples = 0x0;
  
  typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>
    (*(((mode == NFFT_CONV_C2NC) ? in : out )->get_dimensions())); 
  bool oversampled_image = (image_dims==matrix_size_os); 
  
  if( !oversampled_image ){
    throw std::runtime_error("Error: cuNFFT_plan::convolve: ERROR: oversampled image not provided as input.");
  }

  vector<size_t> vec_dims = to_std_vector(matrix_size_os); 
  {
    cuNDArray<complext<REAL> > *image = ((mode == NFFT_CONV_C2NC) ? in : out );
    for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
      vec_dims.push_back(image->get_size(d));
  }

  switch(mode){

  case NFFT_CONV_C2NC:
  	convolve_NFFT_C2NC( in_int, out_int, accumulate );
  	if( dcw_int ) *out_int *= *dcw_int;
    break;
    
  case NFFT_CONV_NC2C:

    // Density compensation
    if( dcw_int ){
      working_samples = new cuNDArray<complext<REAL> >(*in_int);
      *working_samples *= *dcw_int;
    }
    else{
      working_samples = in_int;
    }
    
    _convolve_NFFT_NC2C<REAL,D,ATOMICS>::apply( this, working_samples, out_int, accumulate );
    
    if( dcw_int ){
      delete working_samples; working_samples = 0x0;
    }    
    break;

  default:
    throw std::runtime_error( "Error: cuNFFT_plan::convolve: unknown mode.");
  }

  restore<complext<REAL>, complext<REAL>, REAL>
    (old_device, out, out, out_int, in, in_int, dcw, dcw_int );
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::fft(cuNDArray<complext<REAL> > *data, NFFT_fft_mode mode, bool do_scale )
{
  cuNDArray<complext<REAL> > *data_int = 0x0;
  int old_device;
  
  prepare<complext<REAL>,dummy,dummy>( device, &old_device, data, &data_int );
  
  typename uint64d<D>::Type _dims_to_transform = counting_vec<size_t,D>();
  vector<size_t> dims_to_transform = to_std_vector( _dims_to_transform );
  
  if( mode == NFFT_FORWARDS ){
    cuNDFFT<REAL>::instance()->fft( data_int, &dims_to_transform, do_scale );
  }
  else{
    cuNDFFT<REAL>::instance()->ifft( data_int, &dims_to_transform, do_scale );
  }

  restore<complext<REAL> ,dummy,dummy>(old_device, data, data, data_int);
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::deapodize( cuNDArray<complext<REAL> > *image, bool fourier_domain)
{
  unsigned char components;
  components = _NFFT_FFT;
  check_consistency( 0x0, image, 0x0, components );

  cuNDArray<complext<REAL> > *image_int = 0x0;
  int old_device;
  
  prepare<complext<REAL>,dummy,dummy>(device, &old_device, image, &image_int );

  typename uint64d<D>::Type image_dims = from_std_vector<size_t, D>(*image->get_dimensions()); 
  bool oversampled_image = (image_dims==matrix_size_os); 
  
  if( !oversampled_image ){
    throw std::runtime_error( "Error: cuNFFT_plan::deapodize: ERROR: oversampled image not provided as input.");
  }
  if (fourier_domain){
  	if (!deapodization_filterFFT)
  		deapodization_filterFFT = 	compute_deapodization_filter(true);
  	*image_int *= *deapodization_filterFFT;
  } else {
  	if (!deapodization_filter)
  		deapodization_filter = compute_deapodization_filter(false);
  	*image_int *= *deapodization_filter;
  }
    
  restore<complext<REAL> ,dummy,dummy>(old_device, image, image, image_int);
}

//
// Private class methods
//

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::check_consistency( cuNDArray<complext<REAL> > *samples, cuNDArray<complext<REAL> > *image,
                                                           cuNDArray<REAL> *weights, unsigned char components )
{

  if( !initialized ){
    throw std::runtime_error( "Error: cuNFFT_plan: Unable to proceed without setup.");
  }
  
  if( (components & _NFFT_CONV_C2NC ) && !preprocessed_C2NC ){
    throw std::runtime_error("Error: cuNFFT_plan: Unable to compute NFFT before preprocessing.");
  }
  
  if( (components & _NFFT_CONV_NC2C ) && !(preprocessed_NC2C || (preprocessed_C2NC && ATOMICS ) ) ){
    throw std::runtime_error("Error: cuNFFT_plan: Unable to compute NFFT before preprocessing.");
  }
  
  if( ((components & _NFFT_CONV_C2NC ) || (components & _NFFT_CONV_NC2C )) && !(image && samples) ){
    throw std::runtime_error("Error: cuNFFT_plan: Unable to process 0x0 input/output.");
  }
  
  if( ((components & _NFFT_FFT) || (components & _NFFT_DEAPODIZATION )) && !image ){
    throw std::runtime_error("Error: cuNFFT_plan: Unable to process 0x0 input.");
  }

  if( image->get_number_of_dimensions() < D ){
    throw std::runtime_error("Error: cuNFFT_plan: Number of image dimensions mismatch the plan.");
  }    

  typename uint64d<D>::Type image_dims = from_std_vector<size_t,D>( *image->get_dimensions() );
  bool oversampled_image = (image_dims==matrix_size_os);
  
  if( !((oversampled_image) ? (image_dims == matrix_size_os) : (image_dims == matrix_size) )){
    throw std::runtime_error("Error: cuNFFT_plan: Image dimensions mismatch.");
  }
  
  if( (components & _NFFT_CONV_C2NC ) || (components & _NFFT_CONV_NC2C )){    
    if( (samples->get_number_of_elements() == 0) || (samples->get_number_of_elements() % (number_of_frames*number_of_samples)) ){
      printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %ld.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\n",samples->get_number_of_elements(), number_of_samples, number_of_frames ); fflush(stdout);
      throw std::runtime_error("Error: cuNFFT_plan: The number of samples is not a multiple of #samples/frame x #frames as requested through preprocessing");
    }
    
    unsigned int num_batches_in_samples_array = samples->get_number_of_elements()/(number_of_frames*number_of_samples);
    unsigned int num_batches_in_image_array = 1;

    for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ ){
      num_batches_in_image_array *= image->get_size(d);
    }
    num_batches_in_image_array /= number_of_frames;

    if( num_batches_in_samples_array != num_batches_in_image_array ){
      printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %ld.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\nLeading to %d batches in the samples array.\nThe number of batches in the image array is %d.\n",samples->get_number_of_elements(), number_of_samples, number_of_frames, num_batches_in_samples_array, num_batches_in_image_array ); fflush(stdout);
      throw std::runtime_error("Error: cuNFFT_plan: Number of batches mismatch between samples and image arrays");
    }
  }
  
  if( components & _NFFT_CONV_NC2C ){
    if( weights ){ 
      if( weights->get_number_of_elements() == 0 ||
          !( weights->get_number_of_elements() == number_of_samples || 
             weights->get_number_of_elements() == number_of_frames*number_of_samples) ){
        printf("\ncuNFFT::check_consistency() failed:\n#elements in the samples array: %ld.\n#samples from preprocessing: %d.\n#frames from preprocessing: %d.\n#weights: %ld.\n",samples->get_number_of_elements(), number_of_samples, number_of_frames, weights->get_number_of_elements() ); fflush(stdout);
        throw std::runtime_error("Error: cuNFFT_plan: The number of weights should match #samples/frame x #frames as requested through preprocessing");
      }
    }
  }  
}

template<class REAL, unsigned int D, bool ATOMICS> 
void Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::barebones()
{	
  // These are the fundamental booleans checked before accessing the various member pointers
  initialized = preprocessed_C2NC = preprocessed_NC2C = false;

  // Clear matrix sizes
  clear(matrix_size);
  clear(matrix_size_os);

  // Clear pointers
  trajectory_positions = 0x0;
  tuples_last = bucket_begin = bucket_end = 0x0;

  // and specify the device
  if (cudaGetDevice(&device) != cudaSuccess) {
    throw cuda_error("Error: cuNFFT_plan::barebones:: unable to get device no");
  }
}

template<class REAL, unsigned int D, bool ATOMICS> 
void Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::wipe( NFFT_wipe_mode mode )
{
  // Get current Cuda device
  int old_device;
  if( cudaGetDevice(&old_device) != cudaSuccess ) {
    throw cuda_error("Error: cuNFFT_plan::wipe: unable to get device no");
  }

  if( device != old_device && cudaSetDevice(device) != cudaSuccess) {
    throw cuda_error("Error: cuNFFT_plan::wipe: unable to set device no");
  }

  if( mode==NFFT_WIPE_ALL && initialized ){
    deapodization_filter.reset();
    initialized = false;
  }
    
  if( preprocessed_NC2C ){
    if( tuples_last )  delete tuples_last;
    if( bucket_begin ) delete bucket_begin;
    if( bucket_end )   delete bucket_end;
  }
  
  if( preprocessed_C2NC || preprocessed_NC2C ){
    delete trajectory_positions;
    preprocessed_C2NC = preprocessed_NC2C = false;
  }

  if( device != old_device && cudaSetDevice(old_device) != cudaSuccess) {
    throw cuda_error("Error: cuNFFT_plan::wipe: unable to restore device no");
  }
}

template<class REAL, unsigned int D, bool ATOMICS> 
void Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::compute_beta()
{	
  // Compute Kaiser-Bessel beta paramter according to the formula provided in 
  // Beatty et. al. IEEE TMI 2005;24(6):799-808.
  for( unsigned int d=0; d<D; d++ )
    beta[d] = (M_PI*std::sqrt((W*W)/(alpha[d]*alpha[d])*(alpha[d]-REAL(0.5))*(alpha[d]-REAL(0.5))-REAL(0.8))); 
}

//
// Grid fictitious trajectory with a single sample at the origin
//

template<class REAL, unsigned int D> __global__ void
compute_deapodization_filter_kernel( typename uintd<D>::Type matrix_size_os, typename reald<REAL,D>::Type matrix_size_os_real, 
                                     REAL W, REAL half_W, REAL one_over_W, 
                                     typename reald<REAL,D>::Type beta, complext<REAL> * __restrict__ image_os )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int num_elements = prod(matrix_size_os);

  if( idx <num_elements ){

    // Compute weight from Kaiser-Bessel filter
    const typename uintd<D>::Type cell_pos = idx_to_co<D>(idx, matrix_size_os);

    // Sample position ("origin")
    const vector_td<REAL,D> sample_pos = REAL(0.5)*matrix_size_os_real;

    // Calculate the distance between the cell and the sample
    vector_td<REAL,D> cell_pos_real = vector_td<REAL,D>(cell_pos);
    const typename reald<REAL,D>::Type delta = abs(sample_pos-cell_pos_real);

    // Compute convolution weight. 
    REAL weight; 
    REAL zero = REAL(0);
    vector_td<REAL,D> half_W_vec( half_W );

    if( weak_greater( delta, half_W_vec ) )
      weight = zero;
    else{ 
      weight = KaiserBessel<REAL>( delta, matrix_size_os_real, one_over_W, beta );
      //if( !isfinite(weight) )
      //weight = zero;
    }
    
    // Output weight
    complext<REAL>  result;
    result.vec[0] = weight; 
    result.vec[1] = zero;
    image_os[idx] = result;
  }
}

//
// Function to calculate the deapodization filter
//

template<class REAL, unsigned int D, bool ATOMICS> boost::shared_ptr<cuNDArray<complext<REAL> > >
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::compute_deapodization_filter( bool FFTed)
{
  std::vector<size_t> tmp_vec_os = to_std_vector(matrix_size_os);

 boost::shared_ptr< cuNDArray<complext<REAL> > > filter( new cuNDArray<complext<REAL> >(tmp_vec_os));
  vector_td<REAL,D> matrix_size_os_real = vector_td<REAL,D>(matrix_size_os);
  
  // Find dimensions of grid/blocks.
  dim3 dimBlock( 256 );
  dim3 dimGrid( (prod(matrix_size_os)+dimBlock.x-1)/dimBlock.x );

  // Invoke kernel
  compute_deapodization_filter_kernel<REAL,D><<<dimGrid, dimBlock>>> 
    ( vector_td<unsigned int,D>(matrix_size_os), matrix_size_os_real, W, REAL(0.5)*W, REAL(1)/W, beta, filter->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();
  
  // FFT
  if (FFTed)
  	fft( filter.get(), NFFT_FORWARDS, false );
  else
  	fft( filter.get(), NFFT_BACKWARDS, false );
  // Reciprocal
  reciprocal_inplace(filter.get());
  return filter;
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::compute_NFFT_C2NC( cuNDArray<complext<REAL> > *image, cuNDArray<complext<REAL> > *samples )
{
  // private method - no consistency check. We trust in ourselves.

  // Deapodization
  deapodize( image );
    
  // FFT
  fft( image, NFFT_FORWARDS );

  // Convolution
  convolve( image, samples, 0x0, NFFT_CONV_C2NC );
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::compute_NFFTH_NC2C( cuNDArray<complext<REAL> > *samples, cuNDArray<complext<REAL> > *image )
{
  // private method - no consistency check. We trust in ourselves.

  // Convolution
  convolve( samples, image, 0x0, NFFT_CONV_NC2C );

  // FFT
  fft( image, NFFT_BACKWARDS );
  
  // Deapodization  
  deapodize( image );
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::compute_NFFTH_C2NC( cuNDArray<complext<REAL> > *image, cuNDArray<complext<REAL> > *samples )
{
  // private method - no consistency check. We trust in ourselves.

  // Deapodization
  deapodize( image, true );
 
  // FFT
  fft( image, NFFT_BACKWARDS );

  // Convolution
  convolve( image, samples, 0x0, NFFT_CONV_C2NC );
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::compute_NFFT_NC2C( cuNDArray<complext<REAL> > *samples, cuNDArray<complext<REAL> > *image )
{
  // private method - no consistency check. We trust in ourselves.

  // Convolution
  convolve( samples, image, 0x0, NFFT_CONV_NC2C );
  
  // FFT
  fft( image, NFFT_FORWARDS );
  
  // Deapodization
  deapodize( image, true );
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::convolve_NFFT_C2NC( cuNDArray<complext<REAL> > *image, cuNDArray<complext<REAL> > *samples, bool accumulate )
{
  // private method - no consistency check. We trust in ourselves.
  
  unsigned int num_batches = 1;
  for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
    num_batches *= image->get_size(d);
  num_batches /= number_of_frames;
  
  /*
    Setup grid and threads
  */

  size_t threads_per_block;
  unsigned int max_coils;
    
  threads_per_block = NFFT_THREADS_PER_KERNEL;
  
  if( cudaDeviceManager::Instance()->major_version(device) == 1 ){
    max_coils = NFFT_MAX_COILS_COMPUTE_1x;
  }
  else{
    max_coils = NFFT_MAX_COILS_COMPUTE_2x;
  }
  
  // We can (only) convolve max_coils batches per run due to shared memory issues. 
  unsigned int domain_size_coils_desired = num_batches;
  unsigned int num_repetitions = domain_size_coils_desired/max_coils + 
    ( ((domain_size_coils_desired%max_coils)==0) ? 0 : 1 );
  unsigned int domain_size_coils = (num_repetitions==1) ? domain_size_coils_desired : max_coils;
  unsigned int domain_size_coils_tail = (num_repetitions==1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions-1)*domain_size_coils;

  // Block and Grid dimensions
  dim3 dimBlock( (unsigned int)threads_per_block );
  dim3 dimGrid( (number_of_samples+dimBlock.x-1)/dimBlock.x, number_of_frames );

  // Calculate how much shared memory to use per thread
  size_t bytes_per_thread = domain_size_coils * sizeof( vector_td<REAL,D> );
  size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof( vector_td<REAL,D> );

  unsigned int double_warp_size_power=0;
  unsigned int __tmp = cudaDeviceManager::Instance()->warp_size(device)<<1;
  while(__tmp!=1){
    __tmp>>=1;
    double_warp_size_power++;
  }
  
  vector_td<REAL,D> matrix_size_os_real = vector_td<REAL,D>( matrix_size_os );

  /*
    Invoke kernel
  */

  for( unsigned int repetition = 0; repetition<num_repetitions; repetition++ ){
    NFFT_convolve_kernel<REAL,D>
      <<<dimGrid, dimBlock, ((repetition==num_repetitions-1) ? dimBlock.x*bytes_per_thread_tail : dimBlock.x*bytes_per_thread)>>>
      ( alpha, beta, W, vector_td<unsigned int,D>(matrix_size_os), vector_td<unsigned int,D>(matrix_size_wrap), number_of_samples,
        (repetition==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
        raw_pointer_cast(&(*trajectory_positions)[0]), 
        image->get_data_ptr()+repetition*prod(matrix_size_os)*number_of_frames*domain_size_coils,
        samples->get_data_ptr()+repetition*number_of_samples*number_of_frames*domain_size_coils, 
        double_warp_size_power, REAL(0.5)*W, REAL(1)/(W), accumulate, matrix_size_os_real );

    CHECK_FOR_CUDA_ERROR();    
  }
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::convolve_NFFT_NC2C( cuNDArray<complext<REAL> > *image, cuNDArray<complext<REAL> > *samples, bool accumulate )
{
  _convolve_NFFT_NC2C<REAL,D,ATOMICS>::apply( this, image, samples, accumulate );
}

template<unsigned int D> struct
_convolve_NFFT_NC2C<float,D,true>{ // True: use atomic operations variant
  static bool apply( cuNFFT_plan<float,D,true> *plan, 
                     cuNDArray<complext<float> > *samples, 
                     cuNDArray<complext<float> > *image, 
                     bool accumulate )
  {   
    //
    // Bring in some variables from the plan
    
    unsigned int device = plan->device;
    unsigned int number_of_frames = plan->number_of_frames;
    unsigned int number_of_samples = plan->number_of_samples;
    typename uint64d<D>::Type matrix_size_os = plan->matrix_size_os;
    typename uint64d<D>::Type matrix_size_wrap = plan->matrix_size_wrap;
    typename reald<float,D>::Type alpha = plan->alpha;
    typename reald<float,D>::Type beta = plan->beta;
    float W = plan->W;
    thrust::device_vector< typename reald<float,D>::Type > *trajectory_positions = plan->trajectory_positions;    

    //
    // Atomic operations are only supported in compute model 2.0 and up
    //

    if( cudaDeviceManager::Instance()->major_version(device) == 1 ){
      throw cuda_error("Error: Atomic NC2C NFFT only supported on device with compute model 2.0 or higher");
    }
    
    // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
    if( !((cudaDeviceManager::Instance()->warp_size(device) & (cudaDeviceManager::Instance()->warp_size(device)-1)) == 0 ) ){
      throw cuda_error("cuNFFT: unsupported hardware (warpSize is not a power of two)");
    }
    
    unsigned int num_batches = 1;
    for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
      num_batches *= image->get_size(d);
    num_batches /= number_of_frames;
    
    //
    //  Setup grid and threads
    //
    
    size_t threads_per_block;
    unsigned int max_coils;
    
    threads_per_block = NFFT_THREADS_PER_KERNEL;
    max_coils = NFFT_MAX_COILS_COMPUTE_2x;
    
    // We can (only) convolve domain_size_coils batches per run due to shared memory issues. 
    unsigned int domain_size_coils_desired = num_batches;
    unsigned int num_repetitions = domain_size_coils_desired/max_coils + 
      ( ((domain_size_coils_desired%max_coils)==0) ? 0 : 1 );
    unsigned int domain_size_coils = (num_repetitions==1) ? domain_size_coils_desired : max_coils;
    unsigned int domain_size_coils_tail = (num_repetitions==1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions-1)*domain_size_coils;
    
    // Block and Grid dimensions
    dim3 dimBlock( (unsigned int)threads_per_block ); 
    dim3 dimGrid( (number_of_samples+dimBlock.x-1)/dimBlock.x, number_of_frames );
    
    // Calculate how much shared memory to use per thread
    size_t bytes_per_thread = domain_size_coils * sizeof( vector_td<float,D> );
    size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof( vector_td<float,D> );
    
    unsigned int double_warp_size_power=0, __tmp = cudaDeviceManager::Instance()->warp_size(device)<<1;
    while(__tmp!=1){
      __tmp>>=1;
      double_warp_size_power++;
    }
    
    vector_td<float,D> matrix_size_os_real = vector_td<float,D>( matrix_size_os );
    
    if( !accumulate ){
      clear(image);
    }
    
    //
    // Invoke kernel
    //
    
    for( unsigned int repetition = 0; repetition<num_repetitions; repetition++ ){
      
      NFFT_H_atomic_convolve_kernel<float,D>
        <<<dimGrid, dimBlock, ((repetition==num_repetitions-1) ? dimBlock.x*bytes_per_thread_tail : dimBlock.x*bytes_per_thread)>>>
        ( alpha, beta, W, vector_td<unsigned int,D>(matrix_size_os), vector_td<unsigned int,D>(matrix_size_wrap), number_of_samples,
          (repetition==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils,
          raw_pointer_cast(&(*trajectory_positions)[0]), 
          samples->get_data_ptr()+repetition*number_of_samples*number_of_frames*domain_size_coils,
          image->get_data_ptr()+repetition*prod(matrix_size_os)*number_of_frames*domain_size_coils,
          double_warp_size_power, float(0.5)*W, float(1)/(W), matrix_size_os_real );
    }
    
    CHECK_FOR_CUDA_ERROR();
   
    return true;
  }
};

template<unsigned int D> struct
_convolve_NFFT_NC2C<double,D,true>{ // True: use atomic operations variant
  // Atomics don't exist for doubles, so this gives a compile error if you actually try to use it.
};

template<class REAL, unsigned int D> struct
_convolve_NFFT_NC2C<REAL,D,false>{ // False: use non-atomic operations variant
  static void apply( cuNFFT_plan<REAL,D,false> *plan,
                     cuNDArray<complext<REAL> > *samples, 
                     cuNDArray<complext<REAL> > *image, 
                     bool accumulate )
  {
    // Bring in some variables from the plan
    
    unsigned int device = plan->device;
    unsigned int number_of_frames = plan->number_of_frames;
    unsigned int number_of_samples = plan->number_of_samples;
    typename uint64d<D>::Type matrix_size_os = plan->matrix_size_os;
    typename uint64d<D>::Type matrix_size_wrap = plan->matrix_size_wrap;
    typename reald<REAL,D>::Type alpha = plan->alpha;
    typename reald<REAL,D>::Type beta = plan->beta;
    REAL W = plan->W;
    thrust::device_vector< typename reald<REAL,D>::Type > *trajectory_positions = plan->trajectory_positions;    
    thrust::device_vector<unsigned int> *tuples_last = plan->tuples_last;
    thrust::device_vector<unsigned int> *bucket_begin = plan->bucket_begin;
    thrust::device_vector<unsigned int> *bucket_end = plan->bucket_end;

    // private method - no consistency check. We trust in ourselves.
    // Check if warp_size is a power of two. We do some modulus tricks in the kernels that depend on this...
    if( !((cudaDeviceManager::Instance()->warp_size(device) & (cudaDeviceManager::Instance()->warp_size(device)-1)) == 0 ) ){
      throw cuda_error("cuNFFT: unsupported hardware (warpSize is not a power of two)");

    }
    unsigned int num_batches = 1;
    for( unsigned int d=D; d<image->get_number_of_dimensions(); d++ )
      num_batches *= image->get_size(d);
    num_batches /= number_of_frames;
    
    //
    // Setup grid and threads
    //
    
    size_t threads_per_block;
    unsigned int max_coils;
    
    threads_per_block = NFFT_THREADS_PER_KERNEL;
    
    if( cudaDeviceManager::Instance()->major_version(device) == 1 ){
      max_coils = NFFT_MAX_COILS_COMPUTE_1x;
    }
    else{
      max_coils = NFFT_MAX_COILS_COMPUTE_2x;
    }
    
    // We can (only) convolve domain_size_coils batches per run due to shared memory issues. 
    unsigned int domain_size_coils_desired = num_batches;
    unsigned int num_repetitions = domain_size_coils_desired/max_coils + 
      ( ((domain_size_coils_desired%max_coils)==0) ? 0 : 1 );
    unsigned int domain_size_coils = (num_repetitions==1) ? domain_size_coils_desired : max_coils;
    unsigned int domain_size_coils_tail = (num_repetitions==1) ? domain_size_coils_desired : domain_size_coils_desired - (num_repetitions-1)*domain_size_coils;
    
    // Block and Grid dimensions
    dim3 dimBlock( (unsigned int)threads_per_block ); 
    dim3 dimGrid( (prod(matrix_size_os+matrix_size_wrap)+dimBlock.x-1)/dimBlock.x, number_of_frames );
    
    // Calculate how much shared memory to use per thread
    size_t bytes_per_thread = domain_size_coils * sizeof( vector_td<REAL,D> );
    size_t bytes_per_thread_tail = domain_size_coils_tail * sizeof( vector_td<REAL,D> );
    
    unsigned int double_warp_size_power=0, __tmp = cudaDeviceManager::Instance()->warp_size(device)<<1;
    while(__tmp!=1){
      __tmp>>=1;
      double_warp_size_power++;
    }
    
    vector_td<REAL,D> matrix_size_os_real = vector_td<REAL,D>( matrix_size_os );
    
    // Define temporary image that includes a wrapping zone
    cuNDArray<complext<REAL> > _tmp;
    
    vector<size_t> vec_dims = to_std_vector(matrix_size_os+matrix_size_wrap); 
    if( number_of_frames > 1 )
      vec_dims.push_back(number_of_frames);
    if( num_batches > 1 ) 
      vec_dims.push_back(num_batches);
    
    _tmp.create(&vec_dims);
    
    //
    // Invoke kernel
    //
    
    for( unsigned int repetition = 0; repetition<num_repetitions; repetition++ ){
      
      NFFT_H_convolve_kernel<REAL,D>
        <<<dimGrid, dimBlock, ((repetition==num_repetitions-1) ? dimBlock.x*bytes_per_thread_tail : dimBlock.x*bytes_per_thread)>>>
        ( alpha, beta, W, vector_td<unsigned int,D>(matrix_size_os+matrix_size_wrap), number_of_samples,
          (repetition==num_repetitions-1) ? domain_size_coils_tail : domain_size_coils, 
          raw_pointer_cast(&(*trajectory_positions)[0]), 
          _tmp.get_data_ptr()+repetition*prod(matrix_size_os+matrix_size_wrap)*number_of_frames*domain_size_coils,
          samples->get_data_ptr()+repetition*number_of_samples*number_of_frames*domain_size_coils, 
          raw_pointer_cast(&(*tuples_last)[0]), raw_pointer_cast(&(*bucket_begin)[0]), raw_pointer_cast(&(*bucket_end)[0]),
          double_warp_size_power, REAL(0.5)*W, REAL(1)/(W), matrix_size_os_real );
    }
    
    CHECK_FOR_CUDA_ERROR();
    
    plan->image_wrap( &_tmp, image, accumulate );
  };
};

// Image wrap kernels

template<class REAL, unsigned int D> __global__ void
image_wrap_kernel( typename uintd<D>::Type matrix_size_os, typename uintd<D>::Type matrix_size_wrap, bool accumulate,
                   const complext<REAL> * __restrict__ in, complext<REAL> * __restrict__ out )
{
  unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int num_elements_per_image_src = prod(matrix_size_os+matrix_size_wrap);
  const unsigned int image_offset_src = blockIdx.y*num_elements_per_image_src;
  
  const typename uintd<D>::Type co = idx_to_co<D>(idx, matrix_size_os);
  const typename uintd<D>::Type half_wrap = matrix_size_wrap>>1;
  
  // Make "boolean" vectors denoting whether wrapping needs to be performed in a given direction (forwards/backwards)
  vector_td<bool,D> B_l = vector_less( co, half_wrap );
  vector_td<bool,D> B_r = vector_greater_equal( co, matrix_size_os-half_wrap );
  
  complext<REAL>  result = in[co_to_idx<D>(co+half_wrap, matrix_size_os+matrix_size_wrap) + image_offset_src];

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
              	stride[i-1] = -1;
                wrap_requests++;
              }
              else
              	stride[i-1] = 0;
            }
            else{ 
              if( B_l.vec[i-1] ){ // Wrapping required 
              	stride[i-1] =1 ;
                wrap_requests++;
              }
              else
              	stride[i-1] = 0;
            }
          }
          else{
            // Do not test for wrapping in dimension 'i-1' (for this combination)
          	stride[i-1] = 0;
            skipped_dims++;
          }
        }
	
        // Now it is time to do the actual wrapping (if needed)
        if( wrap_requests == d ){
          typename intd<D>::Type src_co_int = vector_td<int,D>(co+half_wrap);
          typename intd<D>::Type matrix_size_os_int = vector_td<int,D>(matrix_size_os);
          typename intd<D>::Type co_offset_int = src_co_int + component_wise_mul<int,D>(stride,matrix_size_os_int);
          typename uintd<D>::Type co_offset = vector_td<unsigned int,D>(co_offset_int);
          result += in[co_to_idx<D>(co_offset, matrix_size_os+matrix_size_wrap) + image_offset_src];
          break; // only one stride per combination can contribute (e.g. one edge, one corner)
        } 
      } 
    }
  }
  
  // Output
  const unsigned int image_offset_tgt = blockIdx.y*prod(matrix_size_os);
  if( accumulate ) result += out[idx+image_offset_tgt];
  out[idx+image_offset_tgt] = result;
}

template<class REAL, unsigned int D, bool ATOMICS> void
Gadgetron::cuNFFT_plan<REAL,D,ATOMICS>::image_wrap( cuNDArray<complext<REAL> > *source, cuNDArray<complext<REAL> > *target, bool accumulate )
{
  unsigned int num_batches = 1;
  for( unsigned int d=D; d<source->get_number_of_dimensions(); d++ )
    num_batches *= source->get_size(d);
  num_batches /= number_of_frames;

  // Set dimensions of grid/blocks.
  unsigned int bdim = 256;
  dim3 dimBlock( bdim );
  dim3 dimGrid( prod(matrix_size_os)/bdim, number_of_frames*num_batches );

  // Safety check
  if( (prod(matrix_size_os)%bdim) != 0 ) {
  	std::stringstream ss;
  	ss << "Error: cuNFFT : the number of oversampled image elements must be a multiplum of the block size: " << bdim;
    throw std::runtime_error(ss.str());
  }

  // Invoke kernel
  image_wrap_kernel<REAL,D><<<dimGrid, dimBlock>>>
    ( vector_td<unsigned int,D>(matrix_size_os), vector_td<unsigned int,D>(matrix_size_wrap), accumulate, source->get_data_ptr(), target->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();
}	

//
// Template instantion
//

template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 1, true >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 1, false >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< double, 1, false >;

template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 2, true >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 2, false >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< double, 2, false >;

template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 3, true >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 3, false >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< double, 3, false >;

template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 4, true >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< float, 4, false >;
template class EXPORTGPUNFFT Gadgetron::cuNFFT_plan< double, 4, false >;
