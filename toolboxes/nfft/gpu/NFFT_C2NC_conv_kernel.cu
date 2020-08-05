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

template<class REAL> __inline__ __device__ void
NFFT_output( unsigned int number_of_samples, unsigned int number_of_batches, complext<REAL> * __restrict__ samples,
	     unsigned int double_warp_size_power, unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx, bool accumulate )
{
  
  REAL *shared_mem = (REAL*) _shared_mem;
  
  for( unsigned int batch=0; batch<number_of_batches; batch++ ){
    complext<REAL>sample_value;
    sample_value._real = shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)];
    sample_value._imag = shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)+warpSize];

    unsigned int out_idx = (batch*gridDim.y+blockIdx.y)*number_of_samples + globalThreadId;

    if( accumulate ) sample_value += samples[out_idx];
    samples[out_idx] = sample_value;
  }
}

template<unsigned int D> __inline__ __device__ static void
resolve_wrap( vector_td<int,D> &grid_position, vector_td<unsigned int,D> &matrix_size_os )
{
  vector_td<int,D> zero(0);
  grid_position += vector_less(grid_position, zero)*matrix_size_os;
  grid_position -= vector_greater_equal(grid_position, matrix_size_os)* matrix_size_os;
}

template<class REAL, unsigned int D> __inline__ __device__ void
NFFT_iterate_body( typename reald<REAL,D>::Type alpha, typename reald<REAL,D>::Type beta, REAL W, 
		   vector_td<unsigned int, D> matrix_size_os, unsigned int number_of_batches, const complext<REAL> * __restrict__ image,
		   unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,D> matrix_size_os_real, unsigned int sharedMemFirstSampleIdx,
		   vector_td<REAL,D> sample_position, vector_td<int,D> grid_position )
{

    using namespace thrust::cuda_cub;
      
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
  resolve_wrap<D>( grid_position, matrix_size_os);

  REAL *shared_mem = (REAL*) _shared_mem;

  unsigned int image_idx = ((blockIdx.y)*prod(matrix_size_os) + co_to_idx( vector_td<unsigned int, D>(grid_position), matrix_size_os ))*number_of_batches;

  for( unsigned int batch=0; batch<number_of_batches; batch++ ){
    
    // Read the grid cell value from global memory
    const complext<REAL> grid_value = cub::ThreadLoad<cub::LOAD_LDG>(image+image_idx+batch );
    
    // Add 'weight*grid_value' to the samples in shared memory
    shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)] += (weight*grid_value._real);
    shared_mem[sharedMemFirstSampleIdx+(batch<<double_warp_size_power)+warpSize] += (weight*grid_value._imag);
  }
}

//
// This method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,1>::Type alpha, typename reald<REAL,1>::Type beta, REAL W, 
	      vector_td<unsigned int,1> matrix_size_os, unsigned int number_of_batches, const complext<REAL> * __restrict__ image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,1> matrix_size_os_real, unsigned int sharedMemFirstSampleIdx,
	      vector_td<REAL,1> sample_position, vector_td<int,1> lower_limit, vector_td<int,1> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
    
    const intd<1>::Type grid_position(x);
    
    NFFT_iterate_body<REAL,1>( alpha, beta, W, matrix_size_os, number_of_batches, image, double_warp_size_power, half_W, 
			       one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, grid_position );
  }
}

//
// This method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,2>::Type alpha, typename reald<REAL,2>::Type beta, REAL W, 
	      vector_td<unsigned int,2> matrix_size_os, unsigned int number_of_batches, const complext<REAL> * __restrict__ image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,2> matrix_size_os_real, unsigned int sharedMemFirstSampleIdx,
	      vector_td<REAL,2> sample_position, vector_td<int,2> lower_limit, vector_td<int,2> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
    for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
      
      const intd<2>::Type grid_position(x,y);
      
      NFFT_iterate_body<REAL,2>( alpha, beta, W, matrix_size_os, number_of_batches, image, double_warp_size_power, half_W, 
				 one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, grid_position );
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,3>::Type alpha, typename reald<REAL,3>::Type beta, REAL W, 
	      vector_td<unsigned int,3> matrix_size_os, unsigned int number_of_batches, const complext<REAL> * __restrict__ image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,3> matrix_size_os_real, unsigned int sharedMemFirstSampleIdx,
	      vector_td<REAL,3> sample_position, vector_td<int,3> lower_limit, vector_td<int,3> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++ ){
    for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
      for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
	
	const intd<3>::Type grid_position(x,y,z);
	
	NFFT_iterate_body<REAL,3>( alpha, beta, W, matrix_size_os, number_of_batches, image, double_warp_size_power, half_W, 
				   one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, grid_position );
      }
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class REAL> __inline__ __device__ void
NFFT_iterate( typename reald<REAL,4>::Type alpha, typename reald<REAL,4>::Type beta, REAL W, 
	      vector_td<unsigned int,4> matrix_size_os, unsigned int number_of_batches, const complext<REAL> * __restrict__ image,
	      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,4> matrix_size_os_real, unsigned int sharedMemFirstSampleIdx,
	      vector_td<REAL,4> sample_position, vector_td<int,4> lower_limit, vector_td<int,4> upper_limit )
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int w = lower_limit.vec[3]; w<=upper_limit.vec[3]; w++ ){
    for( int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++ ){
      for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
	for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
	  
	  const intd<4>::Type grid_position(x,y,z,w);
	  
	  NFFT_iterate_body<REAL,4>( alpha, beta, W, matrix_size_os, number_of_batches, image, double_warp_size_power, half_W, 
				     one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, grid_position );
	}
      }
    }
  }
}

template<class REAL, unsigned int D> __inline__ __device__ void
NFFT_convolve( typename reald<REAL,D>::Type alpha, typename reald<REAL,D>::Type beta, REAL W, 
	       vector_td<unsigned int, D> matrix_size_os, vector_td<unsigned int, D> matrix_size_wrap, 
	       unsigned int number_of_samples, unsigned int number_of_batches, const vector_td<REAL,D> * __restrict__ traj_positions, const complext<REAL> * __restrict__ image,
	       unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,D> matrix_size_os_real,
	       unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx )
{
  
  // Sample position to convolve onto
  // Computed in preprocessing, which included a wrap zone. Remove this wrapping.
  const vector_td<REAL,D> half_wrap_real = vector_td<REAL,D>(matrix_size_wrap>>1);
  const vector_td<REAL,D> sample_position = traj_positions[globalThreadId+blockIdx.y*number_of_samples]-half_wrap_real;
  
  // Half the kernel width
  const vector_td<REAL,D> half_W_vec( half_W );
  
  // Limits of the subgrid to consider
  const vector_td<int,D> lower_limit = vector_td<int,D>( ceil(sample_position-half_W_vec));
  const vector_td<int,D> upper_limit = vector_td<int,D>( floor(sample_position+half_W_vec));

  // Accumulate contributions from the grid
  NFFT_iterate<REAL>( alpha, beta, W, matrix_size_os, number_of_batches, image, double_warp_size_power, 
		      half_W, one_over_W, matrix_size_os_real, sharedMemFirstSampleIdx, sample_position, lower_limit, upper_limit );
}

//
// kernel main
//

template<class REAL, unsigned int D> __global__ void
NFFT_convolve_kernel( typename reald<REAL,D>::Type alpha, typename reald<REAL,D>::Type beta, REAL W, 
		      vector_td<unsigned int, D> matrix_size_os, vector_td<unsigned int, D> matrix_size_wrap,
		      unsigned int number_of_samples, unsigned int number_of_batches, 
		      const vector_td<REAL,D> * __restrict__ traj_positions, const complext<REAL> * __restrict__ image,  complext<REAL> * __restrict__ samples,
		      unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, bool accumulate, vector_td<REAL,D> matrix_size_os_real )
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
  const REAL zero = REAL(0);

  // Initialize shared memory
  for( unsigned int i=0; i<num_reals; i++ )
    shared_mem[sharedMemFirstSampleIdx+warpSize*i] = zero;
  
  // Compute NFFT using arbitrary sample trajectories
  NFFT_convolve<REAL,D>( alpha, beta, W, matrix_size_os, matrix_size_wrap, number_of_samples, number_of_batches, 
			 traj_positions, image, double_warp_size_power, half_W, one_over_W, 
			 matrix_size_os_real, globalThreadId, sharedMemFirstSampleIdx );
  
  // Output k-space image to global memory
  NFFT_output<REAL>( number_of_samples, number_of_batches, samples, double_warp_size_power, globalThreadId, sharedMemFirstSampleIdx, accumulate );
}


