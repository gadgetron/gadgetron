#pragma once
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


#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
/**
 * \brief Atomic addition for double-precision floating-point types.
 *
 * Atomic addition is not available by default on devices with compute
 * capability lower than 6.0, so it needs to be defined manually.
 *
 * This implementation is based on the CUDA C Programming Guide.
 * 
 * \param address Address of data.
 * \param val Value to be added.
 * \return double Result.
 */
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));

    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif


template<class T, unsigned int D, template<class, unsigned int> class K>
__inline__ __device__
static void NFFT_iterate_body(vector_td<unsigned int, D> matrix_size_os, 
		   unsigned int number_of_batches, const T * __restrict__ samples,  T * __restrict__ image,
		   unsigned int frame, unsigned int num_frames,
		   unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
       vector_td<realType_t<T>,D> sample_position, vector_td<int,D> grid_position,
       const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{
    // Calculate the distance between current sample and the grid cell.
    vector_td<realType_t<T>,D> grid_position_real = vector_td<realType_t<T>,D>(grid_position);
    const vector_td<realType_t<T>,D> delta = abs(sample_position - grid_position_real);

    // Compute convolution weight.
    const realType_t<T> weight = kernel->get(delta);

    // Safety measure.
    if (!isfinite(weight))
        return;

    // Resolve wrapping of grid position
    resolve_wrap<D>( grid_position, matrix_size_os );

    for( unsigned int batch=0; batch<number_of_batches; batch++ )
    {
        // Read the grid sample value from global memory
        T sample_value = samples[sample_idx_in_batch+batch*num_samples_per_batch];
        
        // Determine the grid cell idx
        unsigned int grid_idx = 
            (batch*num_frames+frame)*prod(matrix_size_os) + co_to_idx( vector_td<unsigned int, D>(grid_position), matrix_size_os );

        // Atomic update.
        if constexpr (is_complex_type_v<T>)
        {
            atomicAdd(&(((realType_t<T>*)image)[(grid_idx<<1)+0]), weight*real(sample_value));
            atomicAdd(&(((realType_t<T>*)image)[(grid_idx<<1)+1]), weight*imag(sample_value));
        }
        else
        {
            atomicAdd(&(((realType_t<T>*)image)[(grid_idx<<1)+0]), sample_value);
        }
    }
}

//
// This method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,1> matrix_size_os, 
	      unsigned int number_of_batches, const T * __restrict__ samples, T * __restrict__ image,
	      unsigned int frame, unsigned int num_frames, 
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
        vector_td<realType_t<T>,1> sample_position,
        vector_td<int,1> lower_limit, vector_td<int,1> upper_limit,
        const ConvolutionKernel<realType_t<T>, 1, K>* kernel)
{
    // Iterate through all grid cells influencing the corresponding sample
    for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ )
    {
        const intd<1>::Type grid_position(x);
        
        NFFT_iterate_body<T, 1>(matrix_size_os, number_of_batches, samples, image,
              frame, num_frames,
              num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position,
              kernel);
    }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,2> matrix_size_os, 
	      unsigned int number_of_batches, const T * __restrict__ samples, T * __restrict__ image,
	      unsigned int frame, unsigned int num_frames, 
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
        vector_td<realType_t<T>,2> sample_position,
        vector_td<int,2> lower_limit, vector_td<int,2> upper_limit,
        const ConvolutionKernel<realType_t<T>, 2, K>* kernel)
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
    for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
      
      const intd<2>::Type grid_position(x,y);
      
      NFFT_iterate_body<T, 2>(matrix_size_os, number_of_batches, samples, image,
				 frame, num_frames,
         num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position,
         kernel);
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,3> matrix_size_os, 
	      unsigned int number_of_batches, const T * __restrict__ samples, T * __restrict__ image,
	      unsigned int frame, unsigned int num_frames, 	      
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
        vector_td<realType_t<T>,3> sample_position,
        vector_td<int,3> lower_limit, vector_td<int,3> upper_limit,
        const ConvolutionKernel<realType_t<T>, 3, K>* kernel)
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++ ){
    for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
      for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
	
	const intd<3>::Type grid_position(x,y,z);
	
	NFFT_iterate_body<T, 3>(matrix_size_os, number_of_batches, samples, image,
				   frame, num_frames,
           num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position,
           kernel);
      }
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,4> matrix_size_os, 
	      unsigned int number_of_batches, const T * __restrict__ samples, T * __restrict image,
        unsigned int frame, unsigned int num_frames, 
	      unsigned int num_samples_per_batch, unsigned int sample_idx_in_batch, 
        vector_td<realType_t<T>,4> sample_position,
        vector_td<int,4> lower_limit, vector_td<int,4> upper_limit,
        const ConvolutionKernel<realType_t<T>, 4, K>* kernel)
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int w = lower_limit.vec[3]; w<=upper_limit.vec[3]; w++ ){
    for( int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++ ){
      for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
	for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
	  
	  const intd<4>::Type grid_position(x,y,z,w);
	  
	  NFFT_iterate_body<T, 4>(matrix_size_os, number_of_batches, samples, image,
				     frame, num_frames,
             num_samples_per_batch, sample_idx_in_batch, sample_position, grid_position,
             kernel);
	}
      }
    }
  }
}

//
// kernel main
//

template<class T, unsigned int D, template<class, unsigned int> class K>
__global__ void
NFFT_H_atomic_convolve_kernel(vector_td<unsigned int, D> matrix_size_os, vector_td<unsigned int, D> matrix_size_wrap,
			       unsigned int num_samples_per_frame, unsigned int num_batches, 
			       const vector_td<realType_t<T>,D> * __restrict__ traj_positions, const T * __restrict__ samples, T * __restrict__ image,
             const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{
    const unsigned int sample_idx_in_frame = (blockIdx.x * blockDim.x + threadIdx.x);

    // Check if we are within bounds
    if( sample_idx_in_frame >= num_samples_per_frame )
      return;
        
    const unsigned int frame = blockIdx.y;
    const unsigned int num_frames = gridDim.y;
    const unsigned int num_samples_per_batch = num_samples_per_frame*num_frames ;
    const unsigned int sample_idx_in_batch = sample_idx_in_frame+frame*num_samples_per_frame;
    
    // Sample position computed in preprocessing includes a wrap zone. Remove this wrapping.
    const vector_td<realType_t<T>,D> half_wrap_real =
        vector_td<realType_t<T>,D>(matrix_size_wrap>>1);
    const vector_td<realType_t<T>,D> sample_position =
        traj_positions[sample_idx_in_batch]-half_wrap_real;
    
    // Half the kernel width
    const vector_td<realType_t<T>,D> radius_vec(kernel->get_radius());
    
    // Limits of the subgrid to consider
    const vector_td<int, D> lower_limit(ceil(sample_position - radius_vec));
    const vector_td<int, D> upper_limit(floor(sample_position + radius_vec));

    // Output to the grid.
    NFFT_iterate<T>(matrix_size_os, num_batches, samples, image,
            frame, num_frames, num_samples_per_batch, sample_idx_in_batch, 
            sample_position, lower_limit, upper_limit, kernel);
}
