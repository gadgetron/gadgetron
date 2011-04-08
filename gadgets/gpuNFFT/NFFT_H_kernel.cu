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
NFFT_H_output( unsigned int number_of_batches, NDTYPE *image, 
	       unsigned int double_warp_size_power, unsigned int number_of_domains, 
	       unsigned int globalThreadId, unsigned int sharedMemFirstCellIdx )
{

  REAL *shared_mem = (REAL*) _shared_mem;
  
  for( unsigned int batch=0; batch<number_of_batches; batch++ ){
    NDTYPE cell_coefficient = make_realComplex( shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)], 
						shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)+warpSize] );
    image[batch*number_of_domains+globalThreadId] = cell_coefficient;
  }
}


template <class UINTd, class REALd, class REAL, class NDTYPE> __inline__ __device__ void
NFFT_H_convolve( REAL alpha, REAL beta, REAL W, 
		 UINTd fixed_dims, unsigned int number_of_samples, unsigned int number_of_batches,
		 REALd *traj_positions, NDTYPE *samples, unsigned int *tuples_last, unsigned int *bucket_begin, unsigned int *bucket_end,
		 unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real, unsigned int globalThreadId, UINTd domainPos, unsigned int sharedMemFirstCellIdx )
{

  REAL *shared_mem = (REAL*) _shared_mem;

  // Cell position as reald
  const REALd cell_pos = uintd_to_reald( domainPos ); 

  // Convolve samples onto the domain (shared memory)
  for( unsigned int i=bucket_begin[globalThreadId]; i<bucket_end[globalThreadId]; i++ )
    {
      // Safety precaution TODO
      unsigned int sampleIdx = tuples_last[i];

      // Safety precaution TODO
      REALd sample_pos = traj_positions[sampleIdx];
      
      // Calculate the distance between the cell and the sample
      REALd delta = abs(sample_pos-cell_pos);
      REALd half_W_vec; real_to_reald(half_W, half_W, half_W_vec );
  
      // Check if sample will contribute
      if( weak_greater(delta, half_W_vec ))
	continue;
      
      // Compute convolution weights
      float weight = KaiserBessel( delta, matrix_size_os_real, one_over_W, beta, fixed_dims );
      
      // Safety measure. We have occationally observed a NaN from the KaiserBessel computation
      if( !isfinite(weight) )
	continue;
      
      // Apply Kaiser-Bessel filter to input images
      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
	
	NDTYPE sample_val = samples[sampleIdx+batch*number_of_samples];
	
	// Apply filter to shared memory domain. 
	shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)] += weight*sample_val.x;
	shared_mem[sharedMemFirstCellIdx+(batch<<double_warp_size_power)+warpSize] += weight*sample_val.y;
      }
    }
}

//
// kernel main
//

template <class UINTd, class REALd, class REAL, class NDTYPE> __global__ void
NFFT_H_convolve_kernel( REAL alpha, REAL beta, REAL W,
			UINTd domain_count_grid, UINTd fixed_dims, unsigned int number_of_samples, unsigned int number_of_batches,
			REALd *traj_positions, NDTYPE *image, NDTYPE *samples, unsigned int *tuples_last, unsigned int *bucket_begin, unsigned int *bucket_end,
			unsigned int double_warp_size_power, REAL half_W, REAL one_over_W, REALd matrix_size_os_real )
{
  
  // Global thread index
  const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;

  // Number of domains
  const unsigned int number_of_domains = prod(domain_count_grid);

  // Check if we are within bounds
  if( index >= number_of_domains )
    return;
  
  // Mapped global thread index
  const unsigned int domainIdx = index;// TODO: domainsMap[number_of_domains-index-1].y;

  // Compute global domain position
  const UINTd domainPos = idx_to_co( domainIdx, domain_count_grid );
	
  // Number of cells
  const unsigned int num_reals = number_of_batches<<1;

  // All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size to limit bank conflicts
  const unsigned int scatterSharedMemStart = (threadIdx.x/warpSize)*warpSize;
  const unsigned int scatterSharedMemStartOffset = threadIdx.x&(warpSize-1); // a faster way of saying (threadIdx.x%warpSize) 
  const unsigned int sharedMemFirstCellIdx = scatterSharedMemStart*num_reals + scatterSharedMemStartOffset;

  REAL *shared_mem = (REAL*) _shared_mem;
  REAL zero; get_zero(zero);

  // Initialize shared memory
  for( unsigned int i=0; i<num_reals; i++ )
    shared_mem[sharedMemFirstCellIdx+warpSize*i] = zero;
  
  // Compute NFFT using arbitrary sample trajectories.
  NFFT_H_convolve<UINTd, REALd, REAL, NDTYPE>
    ( alpha, beta, W, fixed_dims, number_of_samples, number_of_batches,
      traj_positions, samples, tuples_last, bucket_begin, bucket_end,
      double_warp_size_power, half_W, one_over_W, matrix_size_os_real, index, domainPos, sharedMemFirstCellIdx );
  
  // Output k-space image to global memory
  NFFT_H_output<REAL,NDTYPE>( number_of_batches, image, double_warp_size_power, number_of_domains, index, sharedMemFirstCellIdx );
}
