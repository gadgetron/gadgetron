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

//
// There is no header file accompanying this kernel, so it makes most sense to read the code/file from the end and upwards
//

//
// Transfer result from shared memory to global memory.
//

template<class REAL, class NDTYPE> __inline__ __device__ void 
NFFT_output( unsigned int number_of_samples, unsigned int number_of_batches,
	     NDTYPE *samples,
	     unsigned int double_warp_size_power,
	     unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx )
{
  
  REAL *shared_mem = (REAL*) _shared_mem;
  
  for( unsigned int batch=0; batch<number_of_batches; batch++ ){
    NDTYPE sample_value = make_realComplex( shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)], 
					    shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)+warpSize] );
    samples[batch*number_of_samples+globalThreadId] = sample_value;
  }
}

template<class INTd> __inline__ __device__ void
resolve_wrap( INTd &grid_position, INTd matrix_size_os )
{
  INTd zero = make_intd(0,0);
  grid_position += (vec_less(grid_position, zero)*matrix_size_os);
  grid_position -= (vec_greater_equal(grid_position, matrix_size_os)*matrix_size_os);
}

template<class INTd, class UINTd, class REALd, class REAL, class NDTYPE> __inline__ __device__ void
NFFT_iterate_body( REAL alpha, REAL beta, REAL W,
		   UINTd matrix_size_os, UINTd fixed_dims,
		   unsigned int number_of_batches, 
		   NDTYPE *image,
		   unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real,
		   unsigned int sharedMemFirstSampleIdx,
		   REALd sample_position, INTd grid_position )
{
      
  // Calculate the distance between current sample and the grid cell
  REALd delta = abs(sample_position-intd_to_reald(grid_position));
  REALd half_W_vec = real_to_reald<REAL, REALd>( half_W, half_W );

  // If cell too distant from sample then move on to the next cell
  if( weak_greater(delta, half_W_vec ))
    return;

  // Compute convolution weight.
  REAL weight = KaiserBessel( delta, matrix_size_os_real, one_over_W, beta, fixed_dims );

  // Safety measure. We have occationally observed a NaN from the KaiserBessel computation
  if( !isfinite(weight) )
    return;

  // Resolve wrapping of grid position
  resolve_wrap( grid_position, *((INTd*)&matrix_size_os) );

  REAL *shared_mem = (REAL*) _shared_mem;
  
  for( unsigned int batch=0; batch<number_of_batches; batch++ ){
    
    // Read the grid cell value from global memory
    NDTYPE grid_value = image[batch*prod(matrix_size_os)+co_to_idx( *((UINTd*)&grid_position), matrix_size_os )];
    
    // Add 'weight*grid_value' to the samples in shared memory
    shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)] += (weight*grid_value.x);
    shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)+warpSize] += (weight*grid_value.y);
  }
}

//
// This method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
//

template <class INTd, class UINTd, class REALd, class REAL, class NDTYPE> __inline__ __device__ void
NFFT_iterate( REAL alpha, REAL beta, REAL W,
	      uint2 matrix_size_os, uint2 fixed_dims,
	      unsigned int number_of_batches, 
	      NDTYPE *image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real,
	      unsigned int sharedMemFirstSampleIdx,
	      REALd sample_position, INTd lower_limit, INTd upper_limit )
{
  
  // Iterate through all grid cells influencing the corresponding sample
  for( int y = lower_limit.y; y<=upper_limit.y; y++ ){
    for( int x = lower_limit.x; x<=upper_limit.x; x++ ){
      
      INTd grid_position = make_intd(x,y);
      
      NFFT_iterate_body( alpha, beta, W, matrix_size_os, fixed_dims, number_of_batches, image, double_warp_size_power, 
			 half_W, one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, grid_position );
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template < class INTd, class UINTd, class REALd, class REAL, class NDTYPE> __inline__ __device__ void
NFFT_iterate( REAL alpha, REAL beta, REAL W,
	      uint3 matrix_size_os, uint3 fixed_dims,
	      unsigned int number_of_batches, 
	      NDTYPE *image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real, 
	      unsigned int sharedMemFirstSampleIdx,
	      REALd sample_position, INTd lower_limit, INTd upper_limit )
{

  // Iterate through all grid cells influencing the corresponding sample
  for( int z = lower_limit.z; z<=upper_limit.z; z++ ){
    for( int y = lower_limit.y; y<=upper_limit.y; y++ ){
      for( int x = lower_limit.x; x<=upper_limit.x; x++ ){
	
	INTd grid_position = make_intd(x,y,z);
	
	NFFT_iterate_body( alpha, beta, W, matrix_size_os, fixed_dims, number_of_batches, image, double_warp_size_power, 
			   half_W, one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, grid_position );
      }
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template <class INTd, class UINTd, class REALd, class REAL, class NDTYPE> __inline__ __device__ void
NFFT_iterate( REAL alpha, REAL beta, REAL W,
	      uint4 matrix_size_os, uint4 fixed_dims,
	      unsigned int number_of_batches, 
	      NDTYPE *image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real, 
	      unsigned int sharedMemFirstSampleIdx,
	      REALd sample_position, INTd lower_limit, INTd upper_limit )
{

  // Iterate through all grid cells influencing the corresponding sample
  for( int w = lower_limit.w; w<=upper_limit.w; w++ ){
    for( int z = lower_limit.z; z<=upper_limit.z; z++ ){
      for( int y = lower_limit.y; y<=upper_limit.y; y++ ){
	for( int x = lower_limit.x; x<=upper_limit.x; x++ ){
	  
	  INTd grid_position = make_intd(x,y,z,w);
	  
	  NFFT_iterate_body( alpha, beta, W, matrix_size_os, fixed_dims, number_of_batches, image, double_warp_size_power, 
			     half_W, one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, grid_position );
	}
      }
    }
  }
}

template <class INTd, class UINTd, class REALd, class REAL, class NDTYPE> __inline__ __device__ void
NFFT_convolve( REAL alpha, REAL beta, REAL W,
	       UINTd matrix_size_os, UINTd matrix_size_wrap, UINTd fixed_dims,
	       unsigned int number_of_batches, 
	       REALd *traj_positions, NDTYPE *image,
	       unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real, UINTd non_fixed_dims,
	       unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx )
{
  
  // Sample position to convolve onto
  // Computed in preprocessing, which included a wrap zone. Remove this wrapping.
  REALd sample_position = traj_positions[globalThreadId]-uintd_to_reald(matrix_size_wrap>>1);
  
  // Half the kernel width
  REALd half_W_vec = real_to_reald<REAL, REALd>( half_W, half_W );
  
  // Limits of the subgrid to consider
  INTd lower_limit = reald_to_intd(ceil(sample_position-half_W_vec*uintd_to_reald(non_fixed_dims)));
  INTd upper_limit = reald_to_intd(floor(sample_position+half_W_vec*uintd_to_reald(non_fixed_dims)));

  // Accumulate contributions from the grid
  NFFT_iterate<INTd, UINTd, REALd, REAL, NDTYPE>( alpha, beta, W, matrix_size_os, fixed_dims, number_of_batches, image, double_warp_size_power, 
						  half_W, one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, lower_limit, upper_limit );
}

//
// kernel main
//

template <class UINTd, class REALd, class REAL, class NDTYPE> __global__ void
NFFT_convolve_kernel( REAL alpha, REAL beta, REAL W,
		      UINTd matrix_size_os, UINTd matrix_size_wrap, UINTd fixed_dims,
		      unsigned int number_of_samples, unsigned int number_of_batches,
		      REALd *traj_positions, NDTYPE *image, NDTYPE *samples,
		      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real, UINTd non_fixed_dims )
{

  // Global thread number	
  const unsigned int globalThreadId = (blockIdx.x*blockDim.x+threadIdx.x);

  // Check if we are within bounds
  if( globalThreadId >= number_of_samples )
    return;
  
  // Number of reals to compute/output per thread
  const unsigned int num_reals = number_of_batches<<1;
  
  // All shared memory reals corresponding to domain 'threadIdx.x' are located in bank threadIdx.x%warp_size to limit bank conflicts
  const unsigned int scatterSharedMemStart = (threadIdx.x/warpSize)*warpSize;
  const unsigned int scatterSharedMemStartOffset = threadIdx.x&(warpSize-1); // a faster way of saying (threadIdx.x%warpSize) 
  const unsigned int sharedMemFirstSampleIdx = scatterSharedMemStart*num_reals + scatterSharedMemStartOffset;

  REAL *shared_mem = (REAL*) _shared_mem;
  REAL zero = get_zero<REAL>();

  // Initialize shared memory
  for( unsigned int i=0; i<num_reals; i++ )
    shared_mem[sharedMemFirstSampleIdx+warpSize*i] = zero;

  // Compute NFFT using arbitrary sample trajectories.
  switch(sizeof(UINTd)){

  case sizeof(uint2):
    NFFT_convolve<int2, uint2, REALd, REAL, NDTYPE>( alpha, beta, W, matrix_size_os, matrix_size_wrap, fixed_dims, number_of_batches, 
						     traj_positions, image, double_warp_size_power, half_W, one_over_W, 
						     matrix_size_os_real, non_fixed_dims, globalThreadId, sharedMemFirstSampleIdx );
    break;
    /*   
  case sizeof(uint3):
    NFFT_convolve<int3, uint3, REALd, REAL, NDTYPE>( alpha, beta, W, matrix_size_os, matrix_size_wrap, fixed_dims, number_of_batches, 
						     traj_positions, image, double_warp_size_power, half_W, one_over_W, 
						     matrix_size_os_real, non_fixed_dims, globalThreadId, sharedMemFirstSampleIdx );
    break;
    
  case sizeof(uint4):
    NFFT_convolve<int4, uint4, REALd, REAL, NDTYPE>( alpha, beta, W, matrix_size_os, matrix_size_wrap, fixed_dims, number_of_batches, 
						     traj_positions, image, double_warp_size_power, half_W, one_over_W, 
						     matrix_size_os_real, non_fixed_dims, globalThreadId, sharedMemFirstSampleIdx );
    break;
    */
  }
    
  // Output k-space image to global memory
  NFFT_output<REAL, NDTYPE>( number_of_samples, number_of_batches, samples, double_warp_size_power, globalThreadId, sharedMemFirstSampleIdx );
}
