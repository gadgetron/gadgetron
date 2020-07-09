/*
  CUDA implementation of the NFFT.

  -----------

  Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
  T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
  IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

  Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
  T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
  IEEE Transactions on Medical Imaging 2009; 28(12): 1974-1985. 

  Notice:
  This version of the code uses atomic writes and thus differs from the two references above.
*/

//
// There is no header file accompanying this kernel, so it makes most sense to read the code/file from the end and upwards
//

//
// First the implementation of the inner-most loop
// 

template<class REAL, unsigned int D> static __inline__ __device__ void
NFFT_iterate_body( typename reald<REAL,D>::Type alpha, typename reald<REAL,D>::Type beta, 
		   REAL W, vector_td<unsigned int, D> matrix_size_os, 
		   unsigned int number_of_batches, const complext<REAL> * __restrict__ samples,  complext<REAL> * __restrict__ image,
		   unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,D> matrix_size_os_real, 
		   unsigned int frame, unsigned int num_frames,
		   unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
		   vector_td<REAL,D> sample_position, vector_td<int,D> grid_position )
{
  // Calculate the distance between current sample and the grid cell
  vector_td<REAL,D> grid_position_real = vector_td<REAL,D>(grid_position);
  const vector_td<REAL,D> delta = abs(sample_position-grid_position_real);
  const vector_td<REAL,D> half_W_vec(half_W );
  
  // If cell too distant from sample then move on to the next cell
  if( weak_greater( delta, half_W_vec ))
    return;

  // Compute convolution weight.
  const REAL weight = KaiserBessel<REAL>( delta, matrix_size_os_real, one_over_W, beta );
  
  // Safety measure. We have occationally observed a NaN from the KaiserBessel computation
  if( !isfinite(weight) )
    return;

  // Resolve wrapping of grid position
  resolve_wrap<D>( grid_position, matrix_size_os );

  for( unsigned int batch=0; batch<number_of_batches; batch++ ){

    // Read the grid sample value from global memory
    complext<REAL> sample_value = samples[sample_idx_in_batch+batch*num_samples_per_batch];
    
    // Determine the grid cell idx
    unsigned int grid_idx = 
      (batch*num_frames+frame)*prod(matrix_size_os) + co_to_idx( vector_td<unsigned int, D>(grid_position), matrix_size_os );

    // Atomic update of real and imaginary component
    atomicAdd( &(((REAL*)image)[(grid_idx<<1)+0]), weight*real(sample_value) );
    atomicAdd( &(((REAL*)image)[(grid_idx<<1)+1]), weight*imag(sample_value) );
  }
}

//
// This method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,1>::Type alpha, typename reald<REAL,1>::Type beta, 
	      REAL W, vector_td<unsigned int,1> matrix_size_os, 
	      unsigned int number_of_batches, const complext<REAL> * __restrict__ samples, complext<REAL> * __restrict__ image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, 
	      vector_td<REAL,1> matrix_size_os_real, 
	      unsigned int frame, unsigned int num_frames, 
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
	      vector_td<REAL,1> sample_position, vector_td<int,1> lower_limit, vector_td<int,1> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
    
    const intd<1>::Type grid_position(x);
    
    NFFT_iterate_body<REAL,1>( alpha, beta, W, matrix_size_os, number_of_batches, samples, image, double_warp_size_power, 
			       half_W, one_over_W, matrix_size_os_real, frame, num_frames,
			       num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position );
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,2>::Type alpha, typename reald<REAL,2>::Type beta, 
	      REAL W, vector_td<unsigned int,2> matrix_size_os, 
	      unsigned int number_of_batches, const complext<REAL> * __restrict__ samples, complext<REAL> * __restrict__ image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, 
	      vector_td<REAL,2> matrix_size_os_real, 
	      unsigned int frame, unsigned int num_frames, 
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
	      vector_td<REAL,2> sample_position, vector_td<int,2> lower_limit, vector_td<int,2> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
    for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
      
      const intd<2>::Type grid_position(x,y);
      
      NFFT_iterate_body<REAL,2>( alpha, beta, W, matrix_size_os, number_of_batches, samples, image, double_warp_size_power, 
				 half_W, one_over_W, matrix_size_os_real, frame, num_frames,
				 num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position );
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,3>::Type alpha, typename reald<REAL,3>::Type beta, 
	      REAL W, vector_td<unsigned int,3> matrix_size_os, 
	      unsigned int number_of_batches, const complext<REAL> * __restrict__ samples, complext<REAL> * __restrict__ image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, 
	      vector_td<REAL,3> matrix_size_os_real, 
	      unsigned int frame, unsigned int num_frames, 	      
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
	      vector_td<REAL,3> sample_position, vector_td<int,3> lower_limit, vector_td<int,3> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++ ){
    for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
      for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
	
	const intd<3>::Type grid_position(x,y,z);
	
	NFFT_iterate_body<REAL,3>( alpha, beta, W, matrix_size_os, number_of_batches, samples, image, double_warp_size_power, 
				   half_W, one_over_W, matrix_size_os_real, frame, num_frames,
				   num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position );
      }
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,4>::Type alpha, typename reald<REAL,4>::Type beta, 
	      REAL W, vector_td<unsigned int,4> matrix_size_os, 
	      unsigned int number_of_batches, const complext<REAL> * __restrict__ samples, complext<REAL> * __restrict image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W,
	      vector_td<REAL,4> matrix_size_os_real, 
	      unsigned int frame, unsigned int num_frames, 
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
	      vector_td<REAL,4> sample_position, vector_td<int,4> lower_limit, vector_td<int,4> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int w = lower_limit.vec[3]; w<=upper_limit.vec[3]; w++ ){
    for( int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++ ){
      for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
	for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
	  
	  const intd<4>::Type grid_position(x,y,z,w);
	  
	  NFFT_iterate_body<REAL,4>( alpha, beta, W, matrix_size_os, number_of_batches, samples, image, double_warp_size_power, 
				     half_W, one_over_W, matrix_size_os_real, frame, num_frames,
				     num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position );
	}
      }
    }
  }
}

//
// kernel main
//

template<class REAL, unsigned int D> __global__ void
NFFT_H_atomic_convolve_kernel( typename reald<REAL,D>::Type alpha, typename reald<REAL,D>::Type beta, REAL W, 
			       vector_td<unsigned int, D> matrix_size_os, vector_td<unsigned int, D> matrix_size_wrap,
			       unsigned int num_samples_per_frame, unsigned int num_batches, 
			       const vector_td<REAL,D> * __restrict__ traj_positions, const complext<REAL> * __restrict__ samples, complext<REAL> * __restrict__ image,
			       unsigned int double_warp_size_power, REAL half_W, REAL one_over_W,
			       vector_td<REAL,D> matrix_size_os_real )
{
  
  // A runtime check will prevent this kernel from being run for compute models 1.x.
  //
  
#if(__CUDA_ARCH__>=200)
    
  const unsigned int sample_idx_in_frame = (blockIdx.x*blockDim.x+threadIdx.x);

  // Check if we are within bounds
  if( sample_idx_in_frame >= num_samples_per_frame )
    return;
      
  const unsigned int frame = blockIdx.y;
  const unsigned int num_frames = gridDim.y;
  const unsigned int num_samples_per_batch = num_samples_per_frame*num_frames ;
  const unsigned int sample_idx_in_batch = sample_idx_in_frame+frame*num_samples_per_frame;
  
  // Sample position computed in preprocessing includes a wrap zone. Remove this wrapping.
  const vector_td<REAL,D> half_wrap_real = vector_td<REAL,D>(matrix_size_wrap>>1);
  const vector_td<REAL,D> sample_position = traj_positions[sample_idx_in_batch]-half_wrap_real;
  
  // Half the kernel width
  const vector_td<REAL,D> half_W_vec = vector_td<REAL,D>( half_W );
  
  // Limits of the subgrid to consider
  const vector_td<int,D> lower_limit = vector_td<int,D>( ceil(sample_position-half_W_vec));
  const vector_td<int,D> upper_limit = vector_td<int,D>( floor(sample_position+half_W_vec));

  // Output to the grid
  NFFT_iterate<REAL>( alpha, beta, W, matrix_size_os, num_batches, samples, image, double_warp_size_power, 
		      half_W, one_over_W, matrix_size_os_real, 
		      frame, num_frames, num_samples_per_batch, sample_idx_in_batch, 
		      sample_position, lower_limit, upper_limit );
#endif
}
