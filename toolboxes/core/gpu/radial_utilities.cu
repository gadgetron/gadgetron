#include "radial_utilities.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include "real_utilities.h"
#include "real_utilities_device.h"
#include "check_CUDA.h"

#include <math_constants.h>
#include <vector>
#include <iostream>

using namespace std;

namespace Gadgetron{

  template<class REAL, unsigned int GOLDEN_RATIO_ANGULAR_STEP_SIZE> __inline__ __device__ REAL get_angle_step_GR();

  template<> __inline__ __device__ float get_angle_step_GR<float,0>(){ return CUDART_PI_F*(3.0f-::sqrtf(5.0f))*0.5f; }   // GR_SMALLEST
  template<> __inline__ __device__ float get_angle_step_GR<float,1>(){ return CUDART_PI_F/((::sqrtf(5.0f)+1.0f)*0.5f); } // GR_ORIGINAL
  template<> __inline__ __device__ double get_angle_step_GR<double,0>(){ return CUDART_PI*(3.0-::sqrt(5.0))*0.5; }       // GR_SMALLEST
  template<> __inline__ __device__ double get_angle_step_GR<double,1>(){ return CUDART_PI/((::sqrt(5.0)+1.0)*0.5); }     // GR_ORIGINAL

  template<class REAL, unsigned int GOLDEN_RATIO_ANGULAR_STEP_SIZE> __global__ void
  compute_radial_trajectory_golden_ratio_2d_kernel( typename reald<REAL,2>::Type *co, REAL angular_offset )
  {
    const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;              

    const REAL samples_per_profile = (REAL) blockDim.x;
    const REAL bias = samples_per_profile * REAL(0.5);
    const REAL sample_idx_on_profile = (REAL)threadIdx.x;
    const REAL profile = (REAL)blockIdx.x;
    const REAL angle_step = get_angle_step_GR<REAL,GOLDEN_RATIO_ANGULAR_STEP_SIZE>();

    REAL cos_angle, sin_angle;
    gad_sincos<REAL>( (profile+angular_offset)*angle_step+get_pi<REAL>(), &sin_angle, &cos_angle );

    typename reald<REAL,2>::Type sample_pos; 
    sample_pos.vec[0] = (sample_idx_on_profile-bias)*cos_angle/samples_per_profile;
    sample_pos.vec[1] = (sample_idx_on_profile-bias)*sin_angle/samples_per_profile;
  
    co[index] = sample_pos;
  }

  template<class REAL> boost::shared_ptr< cuNDArray< typename reald<REAL,2>::Type > >
  compute_radial_trajectory_golden_ratio_2d( unsigned int num_samples_per_profile, unsigned int num_profiles_per_frame, 
                                             unsigned int num_frames, unsigned int profile_offset, GOLDEN_RATIO_ANGULAR_STEP_SIZE mode )
  {
    typedef typename reald<REAL,2>::Type T;
  
    // Get device properties
    int device; cudaGetDevice( &device );
    cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
    const unsigned int warp_size = deviceProp.warpSize;
  
    if( num_samples_per_profile%warp_size ){
      cout << endl << "compute_radial_trajectory_golden_ratio_2d: #samples/profile is not a multiple of the device's warp size." << endl;
      return boost::shared_ptr< cuNDArray<T> >();
    }

    unsigned int number_of_samples_per_frame = num_samples_per_profile * num_profiles_per_frame;

    // Allocate space for result
    vector<size_t> dims; dims.push_back( number_of_samples_per_frame ); dims.push_back( num_frames );
    boost::shared_ptr< cuNDArray<T> > co( new cuNDArray<T>(&dims) );
  
    if(!co.get()){
      cout << endl << "Error:: compute_radial_trajectory_golden_ratio_2d: memory allocation failed." << endl;
      return boost::shared_ptr< cuNDArray<T> >();
    }
  
    // Set dimensions of grid/blocks.
    dim3 dimBlock( num_samples_per_profile );
    dim3 dimGrid( num_profiles_per_frame*num_frames );
  
    // Invoke kernel (nvcc has been protesting heavily on various other ways to do this...)
    if( mode == GR_SMALLEST )
      compute_radial_trajectory_golden_ratio_2d_kernel<REAL,0><<< dimGrid, dimBlock >>> 
        ( co->get_data_ptr(), (REAL)profile_offset );
    else
      compute_radial_trajectory_golden_ratio_2d_kernel<REAL,1><<< dimGrid, dimBlock >>> 
        ( co->get_data_ptr(), (REAL)profile_offset );
    
    CHECK_FOR_CUDA_ERROR();
  
    return co;
  }
  template<class REAL> __global__ static void
  compute_radial_trajectory_variable_angle_2d_kernel( typename reald<REAL,2>::Type *co,REAL* angles, REAL one_over_num_profiles_per_frame, REAL one_over_num_frames )
  {
    const unsigned int index = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;

    const REAL samples_per_profile = (REAL) blockDim.x;
    const REAL bias = samples_per_profile * REAL(0.5);
    const REAL sample_idx_on_profile = (REAL)threadIdx.x;
    const int frame = blockIdx.y;


    typename reald<REAL,2>::Type sample_pos;
    sample_pos.vec[0] = (sample_idx_on_profile-bias)*cos(angles[frame])/samples_per_profile;
    sample_pos.vec[1] = (sample_idx_on_profile-bias)*sin(angles[frame])/samples_per_profile;

    co[index] = sample_pos;
  }

  template<class REAL> boost::shared_ptr< cuNDArray< typename reald<REAL,2>::Type > >
  compute_radial_trajectory_variable_angle_2d(cuNDArray<REAL>* angles, unsigned int num_samples_per_profile, unsigned int num_profiles_per_frame, unsigned int num_frames, REAL angular_offset )
  {
    typedef typename reald<REAL,2>::Type T;

    // Get device properties
    int device; cudaGetDevice( &device );
    cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
    const unsigned int warp_size = deviceProp.warpSize;

    if( num_samples_per_profile%warp_size ){
      cout << endl << "Error:: compute_radial_trajectory_fixed_angle_2d: #samples/profile is not a multiple of the device's warp size." << endl;
      return boost::shared_ptr< cuNDArray<T> >();
    }

    unsigned int number_of_samples_per_frame = num_samples_per_profile * num_profiles_per_frame;

    // Allocate space for result
    vector<size_t> dims;
    dims.push_back( number_of_samples_per_frame );
    dims.push_back( num_frames );

    boost::shared_ptr< cuNDArray<T> > co( new cuNDArray<T>(&dims) );

    // Set dimensions of grid/blocks.
    dim3 dimBlock( num_samples_per_profile );
    dim3 dimGrid( num_profiles_per_frame, num_frames );

    // Invoke kernel
    compute_radial_trajectory_variable_angle_2d_kernel<REAL><<< dimGrid, dimBlock >>> ( co->get_data_ptr(), angles->get_data_ptr(),REAL(1)/(REAL)num_profiles_per_frame, REAL(1)/(REAL)num_frames);

    CHECK_FOR_CUDA_ERROR();

    return co;
  }

  template<class REAL> __global__ void
  compute_radial_trajectory_fixed_angle_2d_kernel( typename reald<REAL,2>::Type *co, REAL one_over_num_profiles_per_frame, REAL one_over_num_frames, REAL angular_offset )
  {
    const unsigned int index = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;

    const REAL samples_per_profile = (REAL) blockDim.x;
    const REAL bias = samples_per_profile * REAL(0.5);
    const REAL sample_idx_on_profile = (REAL)threadIdx.x;
    const REAL lprofile = (REAL)blockIdx.x;
    const REAL frame = (REAL)blockIdx.y;

    REAL cos_angle, sin_angle;
    gad_sincos<REAL>( (lprofile+frame*one_over_num_frames)*one_over_num_profiles_per_frame*get_pi<REAL>()+angular_offset+get_pi<REAL>(), &sin_angle, &cos_angle );

    typename reald<REAL,2>::Type sample_pos;
    sample_pos.vec[0] = (sample_idx_on_profile-bias)*cos_angle/samples_per_profile;
    sample_pos.vec[1] = (sample_idx_on_profile-bias)*sin_angle/samples_per_profile;

    co[index] = sample_pos;
  }


  template<class REAL> boost::shared_ptr< cuNDArray< typename reald<REAL,2>::Type > > 
  compute_radial_trajectory_fixed_angle_2d( unsigned int num_samples_per_profile, unsigned int num_profiles_per_frame, unsigned int num_frames, REAL angular_offset )
  {
    typedef typename reald<REAL,2>::Type T;
  
    // Get device properties
    int device; cudaGetDevice( &device );
    cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
    const unsigned int warp_size = deviceProp.warpSize;
  
    if( num_samples_per_profile%warp_size ){
      cout << endl << "Error:: compute_radial_trajectory_fixed_angle_2d: #samples/profile is not a multiple of the device's warp size." << endl;
      return boost::shared_ptr< cuNDArray<T> >();
    }

    unsigned int number_of_samples_per_frame = num_samples_per_profile * num_profiles_per_frame;

    // Allocate space for result
    vector<size_t> dims; 
    dims.push_back( number_of_samples_per_frame ); 
    dims.push_back( num_frames );
  
    boost::shared_ptr< cuNDArray<T> > co( new cuNDArray<T>(&dims) );
  
    // Set dimensions of grid/blocks.
    dim3 dimBlock( num_samples_per_profile );
    dim3 dimGrid( num_profiles_per_frame, num_frames );
  
    // Invoke kernel
    compute_radial_trajectory_fixed_angle_2d_kernel<REAL><<< dimGrid, dimBlock >>> ( co->get_data_ptr(), REAL(1)/(REAL)num_profiles_per_frame, REAL(1)/(REAL)num_frames, angular_offset );
  
    CHECK_FOR_CUDA_ERROR();
  
    return co;
  }

  // Find the (eight) neighbors to a given radial sample index

  template<class REAL, unsigned int GOLDEN_RATIO_ANGULAR_STEP_SIZE, bool GR> 
  __inline__ __device__ typename reald<REAL,2>::Type
  compute_radial_neighbors( REAL sample_idx_on_profile, REAL angular_offset, REAL alpha, 
                            REAL one_over_radial_oversampling_factor, REAL one_over_num_profiles,
                            REAL bias, REAL samples_per_profile, REAL profile, REAL num_profiles,
                            typename reald<REAL,2>::Type * __restrict__ p1, typename reald<REAL,2>::Type * __restrict__ p2,
                            typename reald<REAL,2>::Type * __restrict__ p3, typename reald<REAL,2>::Type * __restrict__ p4,
                            typename reald<REAL,2>::Type * __restrict__ p5, typename reald<REAL,2>::Type * __restrict__ p6,
                            typename reald<REAL,2>::Type * __restrict__ p7, typename reald<REAL,2>::Type * __restrict__ p8  )
  {
    // The sample positions (scales) can be either of the _local_ indices 'sample_idx_on_profile' or 'samples_per_projection'-'sample_idx_on_profile'
    // Beware of "skewness" around the origin, i.e. +1 sample one one side
    const REAL ctr_scale       = alpha*((sample_idx_on_profile-bias)*one_over_radial_oversampling_factor);
    const REAL ctr_scale_inv   = alpha*((samples_per_profile-sample_idx_on_profile-bias)*one_over_radial_oversampling_factor);
    const REAL prev_scale      = alpha*((sample_idx_on_profile-bias-1)*one_over_radial_oversampling_factor);
    const REAL prev_scale_inv  = alpha*((samples_per_profile-(sample_idx_on_profile-1)-bias)*one_over_radial_oversampling_factor);
    const REAL next_scale      = alpha*((sample_idx_on_profile-bias+1)*one_over_radial_oversampling_factor);
    const REAL next_scale_inv  = alpha*((samples_per_profile-(sample_idx_on_profile+1)-bias)*one_over_radial_oversampling_factor);
  
    // Unit circle position for current projection
    REAL cos_angle, sin_angle;
  
    switch(GR){
    
    case true: // golden ratio
      {
        const REAL angle_step = get_angle_step_GR<REAL,GOLDEN_RATIO_ANGULAR_STEP_SIZE>();
        gad_sincos<REAL>( (profile+angular_offset)*angle_step, &sin_angle, &cos_angle );
      }
      break;	  
    case false: // fixed angle
      {
        gad_sincos<REAL>( profile*one_over_num_profiles*get_pi<REAL>(), &sin_angle, &cos_angle );	}
      break;
    }
  
    // Find the normal to the current projection direction
    typename reald<REAL,2>::Type normal; normal.vec[0] = -sin_angle; normal.vec[1] = cos_angle;
  
    // The position of the idx itself
    typename reald<REAL,2>::Type sample_pos; sample_pos.vec[0] = ctr_scale*cos_angle; sample_pos.vec[1] = ctr_scale*sin_angle;
  
    // The positions of the previous and next sample
    (*p1).vec[0] = prev_scale*cos_angle; (*p1).vec[1] = prev_scale*sin_angle;
    (*p2).vec[0] = next_scale*cos_angle; (*p2).vec[1] = next_scale*sin_angle;
  
    // Initialize remaining points;
    (*p3).vec[0] = (*p4).vec[0] = (*p5).vec[0] = (*p6).vec[0] = (*p7).vec[0] = (*p8).vec[0] = 
      (*p3).vec[1] = (*p4).vec[1] = (*p5).vec[1] = (*p6).vec[1] = (*p7).vec[1] = (*p8).vec[1] = get_max<REAL>(); // far away...
  
    // Run through all projections to find the closests neighbors
  
    for( unsigned int i=0; i<num_profiles; i++ ){
    
      if( i == profile )
        continue;
    
      // Unit circle position projection 'i'
      switch(GR)
        {
        case true:
          {
            const REAL angle_step = get_angle_step_GR<REAL,GOLDEN_RATIO_ANGULAR_STEP_SIZE>();
            gad_sincos<REAL>( ((REAL)i+angular_offset)*angle_step, &sin_angle, &cos_angle );
          }
          break;

        case false:
          {
            gad_sincos<REAL>( (REAL)i*one_over_num_profiles*get_pi<REAL>(), &sin_angle, &cos_angle );
          }
          break;	
        }

      // Determine sample positions on projection
      typename reald<REAL,2>::Type prev_pos_1;  prev_pos_1.vec[0] = prev_scale*cos_angle;      prev_pos_1.vec[1] = prev_scale*sin_angle;
      typename reald<REAL,2>::Type prev_pos_2;  prev_pos_2.vec[0] = prev_scale_inv*cos_angle;  prev_pos_2.vec[1] = prev_scale_inv*sin_angle;
      typename reald<REAL,2>::Type ctr_pos_1;   ctr_pos_1.vec[0]  = ctr_scale*cos_angle;       ctr_pos_1.vec[1]  = ctr_scale*sin_angle;
      typename reald<REAL,2>::Type ctr_pos_2;   ctr_pos_2.vec[0]  = ctr_scale_inv*cos_angle;   ctr_pos_2.vec[1]  = ctr_scale_inv*sin_angle;
      typename reald<REAL,2>::Type next_pos_1;  next_pos_1.vec[0] = next_scale*cos_angle;      next_pos_1.vec[1] = next_scale*sin_angle;
      typename reald<REAL,2>::Type next_pos_2;  next_pos_2.vec[0] = next_scale_inv*cos_angle;  next_pos_2.vec[1] = next_scale_inv*sin_angle;
    
      // The dot product is used to ensure we find a neighbor on each side
      if( dot<REAL,2>(ctr_pos_1-sample_pos, normal) > REAL(0) ){
    
        if( norm_squared<REAL>(ctr_pos_1-sample_pos) < norm_squared<REAL>(*p4-sample_pos) ){
          *p3 = prev_pos_1;
          *p4 = ctr_pos_1;
          *p5 = next_pos_1;
        }
      }
      else{
     
        if( norm_squared<REAL>(ctr_pos_1-sample_pos) < norm_squared<REAL>(*p7-sample_pos) ){
          *p6 = prev_pos_1;
          *p7 = ctr_pos_1;
          *p8 = next_pos_1;
        }
      }
  
      // The dot product is used to ensure we find a neighbor on each side
      if( dot<REAL,2>(ctr_pos_2-sample_pos, normal) >  REAL(0) ){
  
        if( norm_squared<REAL>(ctr_pos_2-sample_pos) < norm_squared<REAL>(*p4-sample_pos) ){
          *p3 = prev_pos_2;
          *p4 = ctr_pos_2;
          *p5 = next_pos_2;
        }
      }
      else{
      
        if( norm_squared<REAL>(ctr_pos_2-sample_pos) < norm_squared<REAL>(*p7-sample_pos) ){
          *p6 = prev_pos_2;
          *p7 = ctr_pos_2;
          *p8 = next_pos_2;
        }
      }
    }
  
    return sample_pos;
  }

  template<class REAL, unsigned int GOLDEN_RATIO_ANGULAR_STEP_SIZE, bool GR> __global__ void
  compute_radial_dcw_2d_kernel( REAL alpha, REAL one_over_radial_oversampling_factor, REAL one_over_num_profiles, REAL angular_offset, REAL *dcw )
  {
    const REAL samples_per_profile = (REAL) (blockDim.x<<1);
    const REAL sample_idx_on_profile = (REAL)(blockIdx.x*blockDim.x+threadIdx.x);
    const REAL num_profiles = (REAL)gridDim.y;
    const REAL profile = (REAL)blockIdx.y;
    const REAL bias = samples_per_profile*REAL(0.5);

    const unsigned int index = blockIdx.y*samples_per_profile + sample_idx_on_profile;
  
    REAL weight;
  
    if( sample_idx_on_profile == blockDim.x ){

      // Special case - center of profile/k-space
      const REAL radius = (alpha*one_over_radial_oversampling_factor)*REAL(0.5);
      const REAL area = radius*radius*get_pi<REAL>();
      weight = area/num_profiles;
    }
    else{
    
      // General case - all neighbors exist
    
      // Compute sample positions for the current sample and all neighbors
      // The ordering of p1..p8 in the call below follows the edge of the "Voronoi polygon"
    
      typename reald<REAL,2>::Type sample_pos;
      typename reald<REAL,2>::Type p1, p2, p3, p4, p5, p6, p7, p8;
    
      sample_pos = compute_radial_neighbors<REAL,GOLDEN_RATIO_ANGULAR_STEP_SIZE,GR>
        ( sample_idx_on_profile, angular_offset, alpha, 
          one_over_radial_oversampling_factor, one_over_num_profiles, bias, samples_per_profile, profile, num_profiles,
          &p1, &p5, &p2, &p3, &p4, &p8, &p7, &p6 );
    
      // Find midpoints of lines from sample_pos to all other points.
      p1 = REAL(0.5)*(sample_pos+p1); // computing "sample_pos+(p1-sample_pos)/2"
      p2 = REAL(0.5)*(sample_pos+p2);
      p3 = REAL(0.5)*(sample_pos+p3);
      p4 = REAL(0.5)*(sample_pos+p4);
      p5 = REAL(0.5)*(sample_pos+p5);
      p6 = REAL(0.5)*(sample_pos+p6);
      p7 = REAL(0.5)*(sample_pos+p7);
      p8 = REAL(0.5)*(sample_pos+p8);
    
      // The weight is determined by the area of the polygon (http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/)
      weight = REAL(0.5)*
        ((p1.vec[0]*p2.vec[1]-p2.vec[0]*p1.vec[1])+
         (p2.vec[0]*p3.vec[1]-p3.vec[0]*p2.vec[1])+
         (p3.vec[0]*p4.vec[1]-p4.vec[0]*p3.vec[1])+
         (p4.vec[0]*p5.vec[1]-p5.vec[0]*p4.vec[1])+
         (p5.vec[0]*p6.vec[1]-p6.vec[0]*p5.vec[1])+
         (p6.vec[0]*p7.vec[1]-p7.vec[0]*p6.vec[1])+
         (p7.vec[0]*p8.vec[1]-p8.vec[0]*p7.vec[1])+
         (p8.vec[0]*p1.vec[1]-p1.vec[0]*p8.vec[1]));                        
    
      if( weight<REAL(0) ) weight *= -REAL(1);
    }
  
    dcw[index] = weight;
  }

  template<class REAL, unsigned int GOLDEN_RATIO_ANGULAR_STEP_SIZE, bool GR> boost::shared_ptr< cuNDArray<REAL> >
  compute_radial_dcw_2d( unsigned int samples_per_profile, unsigned int num_profiles, 
                         REAL alpha, REAL one_over_radial_oversampling_factor, unsigned int profile_offset = 0 )
  {
    if( num_profiles < 4 ){
      cout << endl << "Error:: compute_radial_dcw_<*>_2d: use at least four profiles" << endl;
      return boost::shared_ptr< cuNDArray<REAL> >();
    }
  
    // Get device properties
    int device; cudaGetDevice( &device );
    cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
    const unsigned int warp_size = deviceProp.warpSize;
  
    if( samples_per_profile%2 ){
      cout << endl << "Error:: compute_radial_dcw_<*>_2d: samples/profile must be even." << endl;
      return boost::shared_ptr< cuNDArray<REAL> >();
    }

    if( samples_per_profile%warp_size ){
      cout << endl << "Error:: compute_radial_dcw_<*>_2d: samples/profile number a multiple of the device's warp size." << endl;
      return boost::shared_ptr< cuNDArray<REAL> >();
    }

    unsigned int number_of_samples = samples_per_profile * num_profiles;
  
    // Allocate space for result
    vector<size_t> dims; dims.push_back( number_of_samples );
    boost::shared_ptr< cuNDArray<REAL> > dcw( new cuNDArray<REAL>(&dims) );
  
    if(!dcw.get()){
      cout << endl << "Error:: compute_radial_dcw_<*>_2d: memory allocation failed." << endl;
      return boost::shared_ptr< cuNDArray<REAL> >();
    }
  
    // Set dimensions of grid/blocks. (division by two due to resource limitations)
    dim3 dimBlock( samples_per_profile>>1 );
    dim3 dimGrid( 2, num_profiles );
  
    // Invoke kernel
    compute_radial_dcw_2d_kernel<REAL,GOLDEN_RATIO_ANGULAR_STEP_SIZE,GR><<< dimGrid, dimBlock >>> 
      ( alpha, one_over_radial_oversampling_factor, REAL(1)/(REAL)num_profiles, (REAL)profile_offset, dcw->get_data_ptr() );
  
    CHECK_FOR_CUDA_ERROR();
  
    return dcw;
  }

  template<class REAL> boost::shared_ptr< cuNDArray<REAL> >
  compute_radial_dcw_golden_ratio_2d( unsigned int samples_per_profile, unsigned int num_profiles, 
                                      REAL alpha, REAL one_over_radial_oversampling_factor, unsigned int profile_offset,
                                      GOLDEN_RATIO_ANGULAR_STEP_SIZE mode)
  {
    if( mode == GR_SMALLEST )
      return compute_radial_dcw_2d<REAL,0,true>
        ( samples_per_profile, num_profiles, alpha, one_over_radial_oversampling_factor, profile_offset );
    else if( mode == GR_ORIGINAL )
      return compute_radial_dcw_2d<REAL,1,true>
        ( samples_per_profile, num_profiles, alpha, one_over_radial_oversampling_factor, profile_offset );
    else
      throw std::runtime_error("\ncompute_radial_dcw_golden_ratio_2d() :: unexpected mode\n");
  }

  template<class REAL> boost::shared_ptr< cuNDArray<REAL> >
  compute_radial_dcw_fixed_angle_2d( unsigned int samples_per_profile, unsigned int num_profiles, 
                                     REAL alpha, REAL one_over_radial_oversampling_factor )
  {
    // The golden ratio template type is ignored when the tailing template argument is false
    return compute_radial_dcw_2d<REAL,GR_ORIGINAL,false>
      ( samples_per_profile, num_profiles, alpha, one_over_radial_oversampling_factor );
  }

  //
  // Instantiation
  //

  template boost::shared_ptr< cuNDArray< typename reald<float,2>::Type > > 
  compute_radial_trajectory_fixed_angle_2d<float>( unsigned int, unsigned int, unsigned int, float );

  template boost::shared_ptr< cuNDArray< typename reald<double,2>::Type > > 
  compute_radial_trajectory_fixed_angle_2d<double>( unsigned int, unsigned int, unsigned int, double );

  template boost::shared_ptr< cuNDArray< typename reald<float,2>::Type > > 
  compute_radial_trajectory_golden_ratio_2d<float>( unsigned int, unsigned int, unsigned int, unsigned int, GOLDEN_RATIO_ANGULAR_STEP_SIZE );

  template boost::shared_ptr< cuNDArray< typename reald<double,2>::Type > > 
  compute_radial_trajectory_golden_ratio_2d<double>( unsigned int, unsigned int, unsigned int, unsigned int, GOLDEN_RATIO_ANGULAR_STEP_SIZE );

  template boost::shared_ptr< cuNDArray<float> >compute_radial_dcw_fixed_angle_2d<float>( unsigned int, unsigned int, float, float);
  template boost::shared_ptr< cuNDArray<double> >compute_radial_dcw_fixed_angle_2d<double>( unsigned int, unsigned int, double, double );

  template boost::shared_ptr< cuNDArray<float> >
  compute_radial_dcw_golden_ratio_2d<float>( unsigned int, unsigned int, float, float, unsigned int, GOLDEN_RATIO_ANGULAR_STEP_SIZE );

  template boost::shared_ptr< cuNDArray<double> >
  compute_radial_dcw_golden_ratio_2d<double>( unsigned int, unsigned int, double, double, unsigned int, GOLDEN_RATIO_ANGULAR_STEP_SIZE );
}
