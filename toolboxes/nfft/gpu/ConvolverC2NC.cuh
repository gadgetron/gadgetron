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
*/

//
// There is no header file accompanying this kernel, so it makes most sense to read the code/file from the end and upwards
//

//
// Transfer result from shared memory to global memory.
//


// Reference to shared memory
extern __shared__ char _shared_mem[];

using namespace Gadgetron;

template<class T>
__inline__ __device__
void NFFT_output(
    unsigned int number_of_samples, unsigned int number_of_batches,
    T* __restrict__ samples, unsigned int warp_size_power,
    unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx,
    bool accumulate)
{
    
    realType_t<T>* shared_mem = (realType_t<T>*) _shared_mem;
    
    for (unsigned int batch=0; batch<number_of_batches; batch++)
    {
        T sample_value;

        if constexpr (is_complex_type_v<T>)
        {
            sample_value._real =
                shared_mem[sharedMemFirstSampleIdx + (batch << warp_size_power)];
            sample_value._imag =
                shared_mem[sharedMemFirstSampleIdx + (batch << warp_size_power) + warpSize];
        }
        else
        {
            sample_value =
                shared_mem[sharedMemFirstSampleIdx + (batch << warp_size_power)];
        }

        unsigned int out_idx =
            (batch * gridDim.y + blockIdx.y) * number_of_samples + globalThreadId;

        if (accumulate)
            sample_value += samples[out_idx];

        samples[out_idx] = sample_value;
    }
}

template<unsigned int D>
__inline__ __device__
static void resolve_wrap(vector_td<int, D>& grid_position,
                         vector_td<unsigned int, D>& matrix_size_os)
{
    vector_td<int, D> zero(0);
    grid_position += vector_less(grid_position, zero)*matrix_size_os;
    grid_position -= vector_greater_equal(grid_position, matrix_size_os)* matrix_size_os;
}

template<class T, unsigned int D, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate_body(vector_td<unsigned int, D> matrix_size_os, unsigned int number_of_batches, const T * __restrict__ image,
		   unsigned int warp_size_power, unsigned int sharedMemFirstSampleIdx,
       vector_td<realType_t<T>,D> sample_position, vector_td<int,D> grid_position,
       const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{

    using namespace thrust::cuda_cub;
      
    // Calculate the distance between current sample and the grid cell
    vector_td<realType_t<T>,D> grid_position_real = vector_td<realType_t<T>,D>(grid_position);
    const vector_td<realType_t<T>,D> delta = abs(sample_position-grid_position_real);

    // Compute convolution weight.
    const realType_t<T> weight = kernel->get(delta);

    // Safety measure.
    if (!isfinite(weight))
        return;

    // Resolve wrapping of grid position.
    resolve_wrap<D>(grid_position, matrix_size_os);

    realType_t<T> *shared_mem = (realType_t<T>*) _shared_mem;

    unsigned int image_idx = ((blockIdx.y)*prod(matrix_size_os) + co_to_idx( vector_td<unsigned int, D>(grid_position), matrix_size_os ))*number_of_batches;

    for (unsigned int batch = 0; batch < number_of_batches; batch++)
    {
        // Read the grid cell value from global memory
        const T grid_value = cub::ThreadLoad<cub::LOAD_LDG>(image+image_idx+batch );
        
        // Add 'weight*grid_value' to the samples in shared memory
        if constexpr (is_complex_type_v<T>)
        {
            shared_mem[sharedMemFirstSampleIdx + (batch << warp_size_power)] += weight * grid_value._real;
            shared_mem[sharedMemFirstSampleIdx + (batch << warp_size_power) + warpSize] += weight * grid_value._imag;
        }
        else
        {
            shared_mem[sharedMemFirstSampleIdx + (batch << warp_size_power)] += weight * grid_value;
        }
    }
}

//
// This method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,1> matrix_size_os, unsigned int number_of_batches, const T* __restrict__ image,
	      unsigned int warp_size_power, unsigned int sharedMemFirstSampleIdx,
        vector_td<realType_t<T>,1> sample_position,
        vector_td<int,1> lower_limit, vector_td<int,1> upper_limit,
        const ConvolutionKernel<realType_t<T>, 1, K>* kernel)
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
    
    const intd<1>::Type grid_position(x);
    
    NFFT_iterate_body<T, 1>(matrix_size_os, number_of_batches, image, warp_size_power,
             sharedMemFirstSampleIdx,
             sample_position, grid_position, kernel);
  }
}

//
// This method is deliberately overloaded in 'UINTd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,2> matrix_size_os, unsigned int number_of_batches, const T* __restrict__ image,
	      unsigned int warp_size_power, unsigned int sharedMemFirstSampleIdx,
        vector_td<realType_t<T>,2> sample_position,
        vector_td<int,2> lower_limit, vector_td<int,2> upper_limit,
        const ConvolutionKernel<realType_t<T>, 2, K>* kernel)
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
    for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
      
      const intd<2>::Type grid_position(x,y);
      
      NFFT_iterate_body<T, 2>(matrix_size_os, number_of_batches, image, warp_size_power,  
         sharedMemFirstSampleIdx,
         sample_position, grid_position, kernel);
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,3> matrix_size_os, unsigned int number_of_batches, const T* __restrict__ image,
	      unsigned int warp_size_power, unsigned int sharedMemFirstSampleIdx,
        vector_td<realType_t<T>,3> sample_position,
        vector_td<int,3> lower_limit, vector_td<int,3> upper_limit,
        const ConvolutionKernel<realType_t<T>, 3, K>* kernel)
{
  // Iterate through all grid cells influencing the corresponding sample
  for( int z = lower_limit.vec[2]; z<=upper_limit.vec[2]; z++ ){
    for( int y = lower_limit.vec[1]; y<=upper_limit.vec[1]; y++ ){
      for( int x = lower_limit.vec[0]; x<=upper_limit.vec[0]; x++ ){
	
	const intd<3>::Type grid_position(x,y,z);
	
	NFFT_iterate_body<T, 3>(matrix_size_os, number_of_batches, image, warp_size_power,  
           sharedMemFirstSampleIdx,
           sample_position, grid_position, kernel);
      }
    }
  }
}

//
// This method is deliberately overloaded in 'd' (rather than templetized) to improve performance of the loop iteration
//

template<class T, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_iterate(vector_td<unsigned int,4> matrix_size_os, unsigned int number_of_batches, const T* __restrict__ image,
	      unsigned int warp_size_power, unsigned int sharedMemFirstSampleIdx,
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
	  
	  NFFT_iterate_body<T, 4>(matrix_size_os, number_of_batches, image, warp_size_power,  
             sharedMemFirstSampleIdx,
             sample_position, grid_position, kernel);
	}
      }
    }
  }
}

template<class T, unsigned int D, template<class, unsigned int> class K>
__inline__ __device__
void NFFT_convolve(vector_td<unsigned int, D> matrix_size_os,
         vector_td<unsigned int, D> matrix_size_wrap, 
	       unsigned int number_of_samples, unsigned int number_of_batches, const vector_td<realType_t<T>,D> * __restrict__ traj_positions, const T* __restrict__ image,
	       unsigned int warp_size_power,
         unsigned int globalThreadId, unsigned int sharedMemFirstSampleIdx,
         const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{
  
    // Sample position to convolve onto
    // Computed in preprocessing, which included a wrap zone. Remove this wrapping.
    const vector_td<realType_t<T>,D> half_wrap_real = vector_td<realType_t<T>,D>(matrix_size_wrap>>1);
    const vector_td<realType_t<T>,D> sample_position = traj_positions[globalThreadId+blockIdx.y*number_of_samples]-half_wrap_real;
    
    // Half the kernel width
    const vector_td<realType_t<T>,D> radius_vec(kernel->get_radius());
    
    // Limits of the subgrid to consider
    const vector_td<int, D> lower_limit(ceil(sample_position - radius_vec));
    const vector_td<int, D> upper_limit(floor(sample_position + radius_vec));

    // Accumulate contributions from the grid
    NFFT_iterate<T>(matrix_size_os, number_of_batches, image, warp_size_power, 
                  sharedMemFirstSampleIdx,
                  sample_position, lower_limit, upper_limit, kernel);
}

//
// kernel main
//

template<class T, unsigned int D, template<class, unsigned int> class K>
__global__
void NFFT_convolve_kernel(vector_td<unsigned int, D> matrix_size_os, vector_td<unsigned int, D> matrix_size_wrap,
		      unsigned int number_of_samples, unsigned int number_of_batches, 
		      const vector_td<realType_t<T>,D>* __restrict__ traj_positions, const T* __restrict__ image,  T* __restrict__ samples,
          unsigned int warp_size_power, bool accumulate,
          const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{
    // Global thread number.
    const unsigned int globalThreadId = (blockIdx.x * blockDim.x + threadIdx.x);

    // Check if we are within bounds.
    if (globalThreadId >= number_of_samples)
        return;
  
    // Number of reals to compute per thread.
    const unsigned int num_reals = is_complex_type_v<T> ?
        number_of_batches << 1 : number_of_batches;

    // For complex numbers, we need twice as many real samples per batch.
    if constexpr (is_complex_type_v<T>)
    {
        warp_size_power += 1;
    }
  
    // All shared memory reals corresponding to domain 'threadIdx.x' are located in bank threadIdx.x%warp_size to limit bank conflicts
    const unsigned int scatterSharedMemStart = (threadIdx.x/warpSize)*warpSize;
    const unsigned int scatterSharedMemStartOffset = threadIdx.x&(warpSize-1); // a faster way of saying (threadIdx.x%warpSize) 
    const unsigned int sharedMemFirstSampleIdx = scatterSharedMemStart*num_reals + scatterSharedMemStartOffset;

    // printf("kernel beta: %f\n", kernel->get_beta()[0]);

    // Initialize shared memory.
    realType_t<T> *shared_mem = (realType_t<T>*) _shared_mem;
    const realType_t<T> zero = realType_t<T>(0);
    for (unsigned int i = 0; i < num_reals; i++)
        shared_mem[sharedMemFirstSampleIdx + warpSize * i] = zero;
  
    // Compute NFFT using arbitrary sample trajectories
    NFFT_convolve<T, D>(matrix_size_os, matrix_size_wrap,
        number_of_samples, number_of_batches, 
        traj_positions, image, warp_size_power,
        globalThreadId, sharedMemFirstSampleIdx, kernel);
    
    // Output k-space image to global memory.
    NFFT_output<T>(number_of_samples, number_of_batches, samples,
        warp_size_power, globalThreadId, sharedMemFirstSampleIdx, accumulate);
}
