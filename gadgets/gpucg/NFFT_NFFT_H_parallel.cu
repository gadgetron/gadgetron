/*

	CUDA implementation of the NFFT_H (gridding).


	PARALLEL mode, i.e. each domain contains one sample point but multiple images.
	-------------


	References:
	-----------

	Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
	T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
	IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

*/

template <class UINTd, class FLOATd> inline __device__ void
NFFT_H_output_parallel( unsigned int globalThreadId, UINTd matrix_size_os, UINTd matrix_size_wrap, unsigned int domain_size_coils, unsigned int double_warp_size_power, unsigned int sharedMemFirstCellIdx, cuFloatComplex *out_data )
{
	const unsigned int num_elements = prod(matrix_size_os+matrix_size_wrap);

	// Output grid cells' complex values to global memory. Each thread writes 'domain_size_coils' complex floats.
	// All shared memory floats corresponding to domain 'threadIdx' is located in bank threadIdx%warp_size, i.e. no bank conflicts.

	for( unsigned int image=0; image<domain_size_coils; image++ ){
		cuFloatComplex cell_coefficient = make_cuFloatComplex( shared_mem[sharedMemFirstCellIdx+(image<<double_warp_size_power)], shared_mem[sharedMemFirstCellIdx+(image<<double_warp_size_power)+my_warpSize] );
		out_data[globalThreadId + image*num_elements] = cell_coefficient;
	}
}

// Radial path
template <class UINTd, class FLOATd, char TYPE, bool CORE> inline __device__ void
NFFT_H_parallel( const UINTd domainPos, const FLOATd matrix_size_os, const UINTd bias, const UINTd bias_os, const unsigned int globalThreadId, const unsigned int domain_size_coils, const float alpha, const float W, const float beta, const unsigned int double_warp_size_power, const float one_over_num_encoding_steps, const unsigned int number_of_samples, const unsigned int number_of_strips, const unsigned int sharedMemFirstCellIdx, const uint2 *stripsMap, const uint2 *strips, bool lastDimFixed )
{
	const float one_over_W = 1.0f/W;

	// Lookup the texture coordinate for the first sample strip for this domain (.x) and number of strips (.y)
	uint2 strip_Index_Count = stripsMap[globalThreadId];

#ifdef BUFFER_SAFETY_CHECK
	if( (strip_Index_Count.y>0) && (((strip_Index_Count.x>=number_of_strips) || (strip_Index_Count.y>(UINT_MAX-strip_Index_Count.x)) || (strip_Index_Count.x+strip_Index_Count.y-1>=number_of_strips ))))
		return;
#endif

	// Convolve each sample strip onto the domain (shared memory)

	// Cell position as floatd
	const FLOATd cell_pos = uintd_to_floatd( domainPos );

	for( unsigned int i=0; i<strip_Index_Count.y; i++ ){

		// Get the current sample strip
		uint2 strip = strips[strip_Index_Count.x+i];

		// Iterate through all samples in current sample strip.
		// strip.x : global index for first sample; strip.y : number of consequtive samples.

#ifdef BUFFER_SAFETY_CHECK
		if( (strip.x>=number_of_samples) || (strip.y>(UINT_MAX-strip.x)) || (strip.x+strip.y-1>=number_of_samples) )
			continue;
#endif

		for( unsigned int sampleIdx = strip.x; sampleIdx < strip.x+strip.y; sampleIdx++ ){

			// Compute the radial sample position			
			FLOATd sample_pos;
			switch(TYPE){
			case 0:
				sample_pos = NFFT_H_get_sample_position<FLOATd,CORE>( sampleIdx );
				break;
			case 1:
			case 2:
				sample_pos = cuda_compute_radial_sample_position<UINTd, FLOATd, TYPE>( sampleIdx, alpha, bias, bias_os );
				break;
			default:
				return;
			}

			// Calculate the distance between the cell and the sample
			const FLOATd delta = fabsf(sample_pos-cell_pos);
		
			// Check if sample will contribute
			if( weak_greater(delta, W/2.0f ))
				continue;

			// Compute convolution weights.
			float weight = KaiserBessel( delta, matrix_size_os, one_over_W, beta, lastDimFixed );

			// Safety measure. We have occationally observed a NaN from the KaiserBessel computation
			if( !isfinite(weight) )
				continue;

			// Apply Kaiser-Bessel filter to input images
			for( unsigned int image=0; image<domain_size_coils; image++ ){

				// Get sample value from global memory
				cuFloatComplex sample_val = tex1Dfetch( TEX_SAMPLE_VALUES_NFFT_H(CORE), sampleIdx+image*number_of_samples );

				// Apply filter to shared memory domain. 
				// 'image<<double_warp_size_power' is equavalent to 'image*double_warp_size', i.e. the shared memory offset for each coil.
				shared_mem[sharedMemFirstCellIdx+(image<<double_warp_size_power)] += weight*sample_val.x;
				shared_mem[sharedMemFirstCellIdx+(image<<double_warp_size_power)+my_warpSize] += weight*sample_val.y;
			}
		}
	} 
}

// Generic path
template <class UINTd, class FLOATd, char TYPE, bool CORE> inline __device__ void
NFFT_H_parallel( const UINTd domainPos, const FLOATd matrix_size_os, const UINTd bias, const UINTd bias_os, const unsigned int globalThreadId, const unsigned int domain_size_coils, const float alpha, const float W, const float beta, const unsigned int double_warp_size_power, const unsigned int number_of_samples, const unsigned int sharedMemFirstCellIdx, FLOATd *traj_positions, unsigned int *tuples_last, unsigned int *bucket_begin, unsigned int *bucket_end, bool lastDimFixed )
{
	const float one_over_W = 1.0f/W;
	const FLOATd cell_pos = uintd_to_floatd( domainPos ); // Cell position as floatd

	// Convolve samples onto the domain (shared memory)
	for( unsigned int i=bucket_begin[globalThreadId]; i<bucket_end[globalThreadId]; i++ )
	{
		unsigned int sampleIdx = tuples_last[i];
		FLOATd sample_pos = traj_positions[sampleIdx];

		// Calculate the distance between the cell and the sample
		const FLOATd delta = fabsf(sample_pos-cell_pos);

		// Check if sample will contribute
		if( weak_greater(delta, W/2.0f ))
			continue;

		// Compute convolution weights.
		float weight = KaiserBessel( delta, matrix_size_os, one_over_W, beta, lastDimFixed );

		// Safety measure. We have occationally observed a NaN from the KaiserBessel computation
		if( !isfinite(weight) )
			continue;

		// Apply Kaiser-Bessel filter to input images
		for( unsigned int image=0; image<domain_size_coils; image++ ){

			cuFloatComplex sample_val = tex1Dfetch( TEX_SAMPLE_VALUES_NFFT_H(CORE), sampleIdx+image*number_of_samples );

			// Apply filter to shared memory domain. 
			// 'image<<double_warp_size_power' is equavalent to 'image*double_warp_size', i.e. the shared memory offset for each coil.
			shared_mem[sharedMemFirstCellIdx+(image<<double_warp_size_power)] += weight*sample_val.x;
			shared_mem[sharedMemFirstCellIdx+(image<<double_warp_size_power)+my_warpSize] += weight*sample_val.y;
		}
	}
}

/*
	This is the 'main' kernel for the parallel NFFT_H
*/

// Parallel kernel. 
// Special case where prod(domain_size_grid) is 1.

#ifndef SELF_SCHEDULING

// Radial path
template <class UINTd, class FLOATd, char TYPE, bool CORE> __global__ void
NFFT_H_convolve_kernel_parallel_radial( UINTd domain_count_grid, UINTd matrix_size_os, UINTd matrix_size_wrap, UINTd bias, UINTd bias_os, float alpha, float W, float beta, unsigned int domain_size_coils, unsigned int double_warp_size_power, float one_over_num_encoding_steps, unsigned int number_of_samples, unsigned int number_of_strips, uint2 *stripsMap, uint2 *strips, uint2 *domainsMap, cuFloatComplex *out_data, bool lastDimFixed )
{
	// Global thread index
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;

	// Number of domains
	const unsigned int number_of_domains = prod(domain_count_grid);

	if( index >= number_of_domains )
		return;

	// Mapped global thread index
	unsigned int domainIdx = domainsMap[number_of_domains-index-1].y;

	// Compute global domain position
	const UINTd domainPos = idx_to_co( domainIdx, domain_count_grid );
	
	// Number of cells
	const unsigned int number_of_cellfloats = domain_size_coils<<1;

	// All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size (i.e. no bank conflicts in warp).
	const unsigned int local_idx = threadIdx.x;
	const unsigned int scatterSharedMemStart = (local_idx/my_warpSize)*my_warpSize;
	const unsigned int scatterSharedMemStartOffset = local_idx&(my_warpSize-1); //local_idx%my_warpSize;
	const unsigned int sharedMemFirstCellIdx = scatterSharedMemStart*number_of_cellfloats + scatterSharedMemStartOffset;

	// Initialize shared memory
	for( unsigned int i=0; i<number_of_cellfloats; i++ )
		shared_mem[sharedMemFirstCellIdx+my_warpSize*i] = 0.0f;

	// Compute NFFT_H using arbitrary distributions (lookup trajectory positions in texture). The result goes to 'grid_cells' in shared memory.
	NFFT_H_parallel<UINTd, FLOATd, TYPE, CORE>
	  ( domainPos, uintd_to_floatd(matrix_size_os), bias, bias_os, domainIdx, domain_size_coils, alpha, W, beta, double_warp_size_power, one_over_num_encoding_steps, number_of_samples, number_of_strips, sharedMemFirstCellIdx, stripsMap, strips, lastDimFixed );

	// Output k-space image to global memory
	NFFT_H_output_parallel<UINTd, FLOATd>( domainIdx, matrix_size_os, matrix_size_wrap, domain_size_coils, double_warp_size_power, sharedMemFirstCellIdx, out_data );
}

// Generic path
template <class UINTd, class FLOATd, char TYPE, bool CORE> __global__ void 
NFFT_H_convolve_kernel_parallel_generic( UINTd domain_count_grid, UINTd matrix_size_os, UINTd matrix_size_wrap, UINTd bias, UINTd bias_os, float alpha, float W, float beta, unsigned int domain_size_coils, unsigned int double_warp_size_power, unsigned int number_of_samples, uint2 *domainsMap, FLOATd *traj_positions, unsigned int *tuples_last, unsigned int *bucket_begin, unsigned int *bucket_end, cuFloatComplex *out_data, bool lastDimFixed )
{
	// Global thread index
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;

	// Number of domains
	const unsigned int number_of_domains = prod(domain_count_grid);

	if( index >= number_of_domains )
		return;

	// Mapped global thread index
	unsigned int domainIdx = index;//domainsMap[number_of_domains-index-1].y;

	// Compute global domain position
	const UINTd domainPos = idx_to_co( domainIdx, domain_count_grid );
	
	// Number of cells
	const unsigned int number_of_cellfloats = domain_size_coils<<1;

	// All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size (i.e. no bank conflicts in warp).
	const unsigned int local_idx = threadIdx.x;
	const unsigned int scatterSharedMemStart = (local_idx/my_warpSize)*my_warpSize;
	const unsigned int scatterSharedMemStartOffset = local_idx&(my_warpSize-1); //local_idx%my_warpSize;
	const unsigned int sharedMemFirstCellIdx = scatterSharedMemStart*number_of_cellfloats + scatterSharedMemStartOffset;

	// Initialize shared memory
	for( unsigned int i=0; i<number_of_cellfloats; i++ )
		shared_mem[sharedMemFirstCellIdx+my_warpSize*i] = 0.0f;

	// Compute NFFT_H using arbitrary distributions (lookup trajectory positions in texture). The result goes to 'grid_cells' in shared memory.
	NFFT_H_parallel<UINTd, FLOATd, TYPE, CORE>
		( domainPos, uintd_to_floatd(matrix_size_os), bias, bias_os, domainIdx, domain_size_coils, alpha, W, beta, double_warp_size_power, number_of_samples, sharedMemFirstCellIdx, traj_positions, tuples_last, bucket_begin, bucket_end, lastDimFixed );

	// Output k-space image to global memory
	NFFT_H_output_parallel<UINTd, FLOATd>( domainIdx, matrix_size_os, matrix_size_wrap, domain_size_coils, double_warp_size_power, sharedMemFirstCellIdx, out_data );
}

#else

/*
	Experimental warp scheduling according to the paper:
	Aila and Laine. Understanding the Efficiency of Ray Traversal on GPUs. High-Performance Graphics 2009.
*/
#include <sm_11_atomic_functions.h>

// How many batches to process before retrieving a new job from the global queue.
#define BATCH_SIZE 1

// Threads per warp
#define B BATCH_SIZE*WARP_SIZE

// Radial path
template <class UINTd, class FLOATd, char TYPE, bool CORE> __global__ void
NFFT_H_convolve_kernel_parallel( UINTd domain_count_grid, UINTd matrix_size_os, UINTd matrix_size_wrap, UINTd bias, UINTd bias_os, float alpha, float W, float beta, unsigned int domain_size_coils, unsigned int double_warp_size_power, float one_over_num_encoding_steps, unsigned int number_of_samples, unsigned int number_of_strips, uint2 *stripsMap, uint2 *strips, uint2 *domainsMap, cuFloatComplex *out_data, bool lastDimFixed )
{
	// variables shared by entire warp, place to shared memory
	__shared__ volatile unsigned int nextThreadArray[BLOCKDIM_Y];
	__shared__ volatile unsigned int threadCountArray[BLOCKDIM_Y];

	volatile unsigned int& localPoolNextThread = nextThreadArray[threadIdx.y];
	volatile unsigned int& localPoolThreadCount = threadCountArray[threadIdx.y];

	// Clear 'threadCountArray'
	if (threadIdx.x==0){
		localPoolThreadCount = 0;
	}

	while (true){

		// get threads from global to local pool
		if (localPoolThreadCount==0 && threadIdx.x==0){
			localPoolNextThread = atomicAdd(&globalPoolNextThread, B);
			localPoolThreadCount = B;
		}

		// get rays from local pool
		int myThreadIndex = localPoolNextThread + threadIdx.x;

		if (myThreadIndex >= globalPoolThreadCount)
			return;

		if (threadIdx.x==0){
			localPoolNextThread += WARP_SIZE;
			localPoolThreadCount -= WARP_SIZE; 
		}

		// Execute kernel (without terminating)

		// Global thread index
		const unsigned int index = myThreadIndex;

		// Number of domains
		const unsigned int number_of_domains = prod(domain_count_grid);

		if( index >= number_of_domains )
			continue;

		// Mapped global thread index
		unsigned int domainIdx = domainsMap[number_of_domains-index-1].y;

		// Compute global domain position
		const UINTd domainPos = idx_to_co( domainIdx, domain_count_grid );

		// Number of cells
		const unsigned int number_of_cellfloats = domain_size_coils<<1;

		// All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size (i.e. no bank conflicts in warp).
		const unsigned int local_idx = threadIdx.x+WARP_SIZE*threadIdx.y;
		const unsigned int scatterSharedMemStart = (local_idx/my_warpSize)*my_warpSize;
		const unsigned int scatterSharedMemStartOffset = local_idx&(my_warpSize-1); //local_idx%my_warpSize;
		const unsigned int sharedMemFirstCellIdx = scatterSharedMemStart*number_of_cellfloats + scatterSharedMemStartOffset;

		// Initialize shared memory
		for( unsigned int i=0; i<number_of_cellfloats; i++ )
			shared_mem[sharedMemFirstCellIdx+my_warpSize*i] = 0.0f;

		// Compute NFFT_H using arbitrary distributions (lookup trajectory positions in texture). The result goes to 'grid_cells' in shared memory.
		NFFT_H_parallel<UINTd, FLOATd, TYPE, CORE>
			( domainPos, uintd_to_floatd(matrix_size_os), bias, bias_os, domainIdx, domain_size_coils, alpha, W, beta, double_warp_size_power, one_over_num_encoding_steps, number_of_samples, number_of_strips, sharedMemFirstCellIdx, stripsMap, strips, lastDimFixed );

		// Output k-space image to global memory
		NFFT_H_output_parallel<UINTd, FLOATd>( domainIdx, matrix_size_os, matrix_size_wrap, domain_size_coils, double_warp_size_power, sharedMemFirstCellIdx, out_data );
	}
}

#endif
