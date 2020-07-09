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
NFFT_H_output( unsigned int number_of_batches, complext<REAL>* __restrict__ image,
	       unsigned int double_warp_size_power, unsigned int number_of_domains, 
	       unsigned int globalThreadId, unsigned int sharedMemFirstCellIdx )
{

  REAL *shared_mem = (REAL*) _shared_mem;
  
  for( unsigned int batch=0; batch<number_of_batches; batch++ ){
    complext<REAL>cell_coefficient;
    cell_coefficient._real = shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)];
    cell_coefficient._imag = shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)+warpSize];
    image[(batch*gridDim.y+blockIdx.y)*number_of_domains+globalThreadId] = cell_coefficient;
  }
}


template<class REAL, unsigned int D> __inline__ __device__ void
NFFT_H_convolve( typename reald<REAL,D>::Type alpha, typename reald<REAL,D>::Type beta, REAL W, 
		 unsigned int number_of_samples, unsigned int number_of_batches, unsigned int number_of_domains,
		 const vector_td<REAL,D> * __restrict__ traj_positions, const  complext<REAL>* __restrict__ samples, const unsigned int * __restrict__ tuples_last,
		 const unsigned int * __restrict__ bucket_begin, const unsigned int * __restrict__ bucket_end,
		 unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, vector_td<REAL,D> matrix_size_os_real, 
		 unsigned int globalThreadId, vector_td<unsigned int,D> domainPos, unsigned int sharedMemFirstCellIdx )
{
    using namespace thrust::cuda_cub;

  REAL *shared_mem = (REAL*) _shared_mem;

  // Cell position as reald
  vector_td<REAL,D> cell_pos = vector_td<REAL,D>( domainPos );
  
  // Convolve samples onto the domain (shared memory)
  const unsigned int frame_offset = blockIdx.y*number_of_domains;
  for( unsigned int i=bucket_begin[globalThreadId+frame_offset]; i<bucket_end[globalThreadId+frame_offset]; i++ )
    {
      // Safety precaution TODO
      unsigned int sampleIdx = tuples_last[i];

      // Safety precaution TODO
      vector_td<REAL,D> sample_pos = traj_positions[sampleIdx];
      
      // Calculate the distance between the cell and the sample
      vector_td<REAL,D> delta = abs(sample_pos-cell_pos);
      vector_td<REAL,D> half_W_vec( half_W );
  
      // Check if sample will contribute
      if( weak_greater(delta, half_W_vec ))
	continue;
      
      // Compute convolution weights
      float weight = KaiserBessel<REAL>( delta, matrix_size_os_real, one_over_W, beta );
      
      // Safety measure
      if( !isfinite(weight) )
      	continue;
      
      // Apply Kaiser-Bessel filter to input images
      for( int batch=0; batch<number_of_batches; batch++ ){
	
//	complext<REAL>sample_val = samples[sampleIdx+batch*gridDim.y*number_of_samples];
       unsigned int idx = sampleIdx*number_of_batches+batch;
	// Apply filter to shared memory domain.

	    complext<REAL> sample_val  = cub::ThreadLoad<cub::LOAD_CS>(samples+idx);
	    shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)] += (weight*sample_val._real);
	    shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)+warpSize] += (weight*sample_val._imag);
      }
    }
}

//
// kernel main
//

template<class REAL, unsigned int D> __global__ void
NFFT_H_convolve_kernel( typename reald<REAL,D>::Type alpha, typename reald<REAL,D>::Type beta, REAL W,
			vector_td<unsigned int,D> domain_count_grid, unsigned int number_of_samples, unsigned int number_of_batches,
			const vector_td<REAL,D> * __restrict__ traj_positions, complext<REAL>* __restrict__ image, const complext<REAL>* __restrict__ samples,
			const unsigned int * __restrict__ tuples_last, const unsigned int * __restrict__ bucket_begin, const unsigned int * __restrict__ bucket_end,
			unsigned int double_warp_size_power,
			REAL half_W, REAL one_over_W, vector_td<REAL,D> matrix_size_os_real )
{
  
  // Global thread index
  const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;

  // Number of domains
  const unsigned int number_of_domains = prod(domain_count_grid);

  // Check if we are within bounds
  if( index >= number_of_domains )
    return;
  
  // Mapped global thread index (actually we don't use a map currently)
  const unsigned int domainIdx = index; 

  // Compute global domain position
  const vector_td<unsigned int,D> domainPos = idx_to_co( domainIdx, domain_count_grid );
	
  // Number of cells
  const unsigned int num_reals = number_of_batches<<1;

  // All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size to limit bank conflicts
  const unsigned int scatterSharedMemStart = (threadIdx.x/warpSize)*warpSize;
  const unsigned int scatterSharedMemStartOffset = threadIdx.x&(warpSize-1); // a faster way of saying (threadIdx.x%warpSize)
  const unsigned int sharedMemFirstCellIdx = scatterSharedMemStart*num_reals + scatterSharedMemStartOffset;

  REAL *shared_mem = (REAL*) _shared_mem;

  // Initialize shared memory
  for( unsigned int i=0; i<num_reals; i++ )
    shared_mem[sharedMemFirstCellIdx+warpSize*i] = REAL(0);
  
  // Compute NFFT using arbitrary sample trajectories.
  NFFT_H_convolve<REAL, D>
    ( alpha, beta, W, number_of_samples, number_of_batches, number_of_domains,
      traj_positions, samples, tuples_last, bucket_begin, bucket_end,
      double_warp_size_power, half_W, one_over_W,  matrix_size_os_real, index, domainPos, sharedMemFirstCellIdx );
  
  // Output k-space image to global memory
  NFFT_H_output<REAL>( number_of_batches, image, double_warp_size_power, number_of_domains, index, sharedMemFirstCellIdx );
}
