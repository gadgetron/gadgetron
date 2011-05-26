#include "radial_utilities.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include "real_utilities.h"
#include "check_CUDA.h"

#include <math_constants.h>
#include <vector>
#include <iostream>

using namespace std;

template<class REAL> __inline__ __device__ REAL get_angle_step_GR();
template<> __inline__ __device__ float get_angle_step_GR(){ return CUDART_PI_F/((sqrtf(5.0f)+1.0f)*0.5f); }
template<> __inline__ __device__ double get_angle_step_GR(){ return CUDART_PI/((sqrt(5.0)+1.0)*0.5); }

template<class REAL> __global__ void
compute_radial_trajectory_golden_ratio_2d_kernel( typename reald<REAL,2>::Type *co, REAL angular_offset )
{
  const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;              

  const REAL samples_per_profile = (REAL) blockDim.x;
  const REAL bias = samples_per_profile * get_half<REAL>();
  const REAL sample_idx_on_profile = (REAL)threadIdx.x;
  const REAL profile = (REAL)blockIdx.x;
  const REAL angle_step = get_angle_step_GR<REAL>();

  REAL cos_angle, sin_angle;
  gad_sincos<REAL>( (profile+angular_offset)*angle_step, &sin_angle, &cos_angle );

  typename reald<REAL,2>::Type sample_pos; 
  sample_pos.vec[0] = (sample_idx_on_profile-bias)*cos_angle/samples_per_profile;
  sample_pos.vec[1] = (sample_idx_on_profile-bias)*sin_angle/samples_per_profile;
  
  co[index] = sample_pos;
}


template<class REAL> auto_ptr< cuNDArray< typename reald<REAL,2>::Type > > 
compute_radial_trajectory_golden_ratio_2d( unsigned int num_samples_per_profile, unsigned int num_profiles_per_frame, unsigned int num_frames, 
					   unsigned int profile_offset )
{
  typedef typename reald<REAL,2>::Type T;

  // Get device properties
  int device; cudaGetDevice( &device );
  cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
  const unsigned int warp_size = deviceProp.warpSize;
  
  if( num_samples_per_profile%warp_size ){
    cout << endl << "compute_radial_trajectory_golden_ratio_2d: #samples/profile number a multiple of the device's warp size." << endl;
    return auto_ptr< cuNDArray<T> >(0x0);
  }

  unsigned int number_of_samples_per_frame = num_samples_per_profile * num_profiles_per_frame;

  // Allocate space for result
  vector<unsigned int> dims; dims.push_back( number_of_samples_per_frame ); dims.push_back( num_frames );
  cuNDArray<T> *co = cuNDArray<T>::allocate(dims);
  
  if(!co){
    cout << endl << "compute_radial_trajectory_golden_ratio_2d: memory allocation failed." << endl;
    return auto_ptr< cuNDArray<T> >(0x0);
  }
  
  // Set dimensions of grid/blocks.
  dim3 dimBlock( num_samples_per_profile );
  dim3 dimGrid( num_profiles_per_frame*num_frames );
  
  // Invoke kernel
  compute_radial_trajectory_golden_ratio_2d_kernel<REAL><<< dimGrid, dimBlock >>> ( co->get_data_ptr(), (REAL)profile_offset );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<T> >(co);  
}

// Find the (eight) neighbors to a given radial sample index

template<class REAL, bool GR> __inline__ __device__ typename reald<REAL,2>::Type
compute_radial_neighbors( REAL sample_idx_on_profile, REAL angular_offset, REAL alpha, 
			  REAL one_over_radial_oversampling_factor, REAL bias, REAL samples_per_profile, REAL profile, REAL num_profiles,
			  typename reald<REAL,2>::Type *p1, typename reald<REAL,2>::Type *p2, typename reald<REAL,2>::Type *p3, typename reald<REAL,2>::Type *p4, 
			  typename reald<REAL,2>::Type *p5, typename reald<REAL,2>::Type *p6, typename reald<REAL,2>::Type *p7, typename reald<REAL,2>::Type *p8  )
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
      const REAL angle_step = get_angle_step_GR<REAL>();
      gad_sincos<REAL>( (profile+angular_offset)*angle_step, &sin_angle, &cos_angle );
    }
    break;
    /*	  
	  case false: // fixed angle
	  {
	  const REAL cur_proj = uint2REAL(profile);
	  const REAL frame = floorf((cur_proj+0.5f)*__one_over_num_profiles);
	  const REAL cur_proj_local = cur_proj - frame*__num_profiles;
	  const REAL rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
	  __sincosf( cur_proj_local*CUDART_PI_F*__one_over_num_profiles+rotation*__interframe_rotation, &sin_angle, &cos_angle );
	  }
	  break;
    */
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
  
  for( unsigned int i=0; i<gridDim.x; i++ ){
    
    if( i == blockIdx.x )
      continue;
    
    // Unit circle position projection 'i'
    switch(GR)
      {
	
      case true:
	{
	  const REAL angle_step = get_angle_step_GR<REAL>();
	  gad_sincos<REAL>( ((REAL)i+angular_offset)*angle_step, &sin_angle, &cos_angle );
	}
	break;
	/*	
      case false:
	{
	  const REAL cur_proj = uint2float(i);
	  const REAL frame = floorf((cur_proj+0.5f)*__one_over_num_profiles);
	  const REAL cur_proj_local = cur_proj - frame*__num_profiles;
	  const REAL rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
	  __sincosf( cur_proj_local*CUDART_PI_F*__one_over_num_profiles+rotation*__interframe_rotation, &sin_angle, &cos_angle );
	}
	break;
	*/
      }
    
    // Determine sample positions on projection
    typename reald<REAL,2>::Type prev_pos_1;  prev_pos_1.vec[0] = prev_scale*cos_angle;      prev_pos_1.vec[1] = prev_scale*sin_angle;
    typename reald<REAL,2>::Type prev_pos_2;  prev_pos_2.vec[0] = prev_scale_inv*cos_angle;  prev_pos_2.vec[1] = prev_scale_inv*sin_angle;
    typename reald<REAL,2>::Type ctr_pos_1;   ctr_pos_1.vec[0]  = ctr_scale*cos_angle;       ctr_pos_1.vec[1]  = ctr_scale*sin_angle;
    typename reald<REAL,2>::Type ctr_pos_2;   ctr_pos_2.vec[0]  = ctr_scale_inv*cos_angle;   ctr_pos_2.vec[1]  = ctr_scale_inv*sin_angle;
    typename reald<REAL,2>::Type next_pos_1;  next_pos_1.vec[0] = next_scale*cos_angle;      next_pos_1.vec[1] = next_scale*sin_angle;
    typename reald<REAL,2>::Type next_pos_2;  next_pos_2.vec[0] = next_scale_inv*cos_angle;  next_pos_2.vec[1] = next_scale_inv*sin_angle;
    
    // The dot product is used to ensure we find a neighbor on each side
    if( dot<REAL,2>(ctr_pos_1-sample_pos, normal) > get_zero<REAL>() ){
    
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
  if( dot<REAL,2>(ctr_pos_2-sample_pos, normal) > get_zero<REAL>() ){
  
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


template<class REAL> __global__ void
compute_radial_dcw_golden_ratio_2d_kernel( REAL alpha, REAL one_over_radial_oversampling_factor, REAL angular_offset, REAL *dcw )
{
  const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;

  const REAL samples_per_profile = (REAL) blockDim.x;
  const REAL num_profiles = (REAL)gridDim.x;
  const REAL bias = samples_per_profile * get_half<REAL>();
  const REAL sample_idx_on_profile = (REAL)threadIdx.x;
  const REAL profile = (REAL)blockIdx.x;
  
  REAL weight;
  
  if( threadIdx.x == (blockDim.x>>1) ){

    // Special case - center of profile/k-space
    const REAL radius = (alpha*one_over_radial_oversampling_factor)*get_half<REAL>();
    const REAL area = radius*radius*get_pi<REAL>();
    weight = area/num_profiles;
  }
  else{
    
    // General case - all neighbors exist
    
    // Compute sample positions for the current sample and all neighbors
    // The ordering of p1..p8 in the call below follows the edge of the "Voronoi polygon"
    
    typename reald<REAL,2>::Type sample_pos;
    typename reald<REAL,2>::Type p1, p2, p3, p4, p5, p6, p7, p8;
    
    sample_pos = compute_radial_neighbors<REAL,true>( sample_idx_on_profile, angular_offset, alpha, 
						      one_over_radial_oversampling_factor, bias, samples_per_profile, profile, num_profiles,
						      &p1, &p5, &p2, &p3, &p4, &p8, &p7, &p6 );
    
    // Find midpoints of lines from sample_pos to all other points.
    p1 = get_half<REAL>()*(sample_pos+p1); // computing "sample_pos+(p1-sample_pos)/2"
    p2 = get_half<REAL>()*(sample_pos+p2);
    p3 = get_half<REAL>()*(sample_pos+p3);
    p4 = get_half<REAL>()*(sample_pos+p4);
    p5 = get_half<REAL>()*(sample_pos+p5);
    p6 = get_half<REAL>()*(sample_pos+p6);
    p7 = get_half<REAL>()*(sample_pos+p7);
    p8 = get_half<REAL>()*(sample_pos+p8);
    
    // The weight is determined by the area of the polygon (http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/)
    weight = get_half<REAL>()*
      ((p1.vec[0]*p2.vec[1]-p2.vec[0]*p1.vec[1])+
       (p2.vec[0]*p3.vec[1]-p3.vec[0]*p2.vec[1])+
       (p3.vec[0]*p4.vec[1]-p4.vec[0]*p3.vec[1])+
       (p4.vec[0]*p5.vec[1]-p5.vec[0]*p4.vec[1])+
       (p5.vec[0]*p6.vec[1]-p6.vec[0]*p5.vec[1])+
       (p6.vec[0]*p7.vec[1]-p7.vec[0]*p6.vec[1])+
       (p7.vec[0]*p8.vec[1]-p8.vec[0]*p7.vec[1])+
       (p8.vec[0]*p1.vec[1]-p1.vec[0]*p8.vec[1]));                        
    
    if( weight<get_zero<REAL>() ) weight *= -get_one<REAL>();
  }
  
  dcw[index] = weight;
}


template<class REAL> auto_ptr< cuNDArray<REAL> >
compute_radial_dcw_golden_ratio_2d( unsigned int samples_per_profile, unsigned int num_profiles, 
				    REAL alpha, REAL one_over_radial_oversampling_factor, unsigned int profile_offset )
{
  if( num_profiles < 4 ){
    cout << endl << "compute_radial_dcw_golden_ratio_2d: use at least four profiles" << endl;
    return auto_ptr< cuNDArray<REAL> >(0x0);
  }
  
  // Get device properties
  int device; cudaGetDevice( &device );
  cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );
  const unsigned int warp_size = deviceProp.warpSize;
  
  if( samples_per_profile%warp_size ){
    cout << endl << "compute_radial_dcw_golden_ratio_2d: samples/profile number a multiple of the device's warp size." << endl;
    return auto_ptr< cuNDArray<REAL> >(0x0);
  }

  unsigned int number_of_samples = samples_per_profile * num_profiles;
  
  // Allocate space for result
  vector<unsigned int> dims; dims.push_back( number_of_samples );
  cuNDArray<REAL> *dcw = cuNDArray<REAL>::allocate(dims);
  
  if(!dcw){
    cout << endl << "compute_radial_dcw_golden_ratio_2d: memory allocation failed." << endl;
    return auto_ptr< cuNDArray<REAL> >(0x0);
  }
  
  // Set dimensions of grid/blocks.
  dim3 dimBlock( samples_per_profile );
  dim3 dimGrid( num_profiles );
  
  // Invoke kernel
  compute_radial_dcw_golden_ratio_2d_kernel<REAL><<< dimGrid, dimBlock >>> ( alpha, one_over_radial_oversampling_factor, (REAL)profile_offset, dcw->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<REAL> >(dcw);  
}


//
// Instantiation
//

template auto_ptr< cuNDArray< typename reald<float,2>::Type > > 
compute_radial_trajectory_golden_ratio_2d<float>( unsigned int, unsigned int, unsigned int, unsigned int );

template auto_ptr< cuNDArray< typename reald<double,2>::Type > > 
compute_radial_trajectory_golden_ratio_2d<double>( unsigned int, unsigned int, unsigned int, unsigned int );

template auto_ptr< cuNDArray<float> >compute_radial_dcw_golden_ratio_2d<float>( unsigned int, unsigned int, float, float, unsigned int );
template auto_ptr< cuNDArray<double> >compute_radial_dcw_golden_ratio_2d<double>( unsigned int, unsigned int, double, double, unsigned int );
