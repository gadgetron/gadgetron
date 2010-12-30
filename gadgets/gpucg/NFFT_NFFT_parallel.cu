/*

	CUDA implementation of the NFFT (inverse gridding).


	PARALLEL mode, i.e. each domain contains one sample point but multiple images.
	-------------


	References:
	-----------

	Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
	T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
	IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

*/


/*
    Block configuration for kernel
	------------------------------
    x,y -> w threads means: shared memory = x, registers = y, max threads = w
	This information is obtained from the .cubin file (x,y) and the Occupancy Calculator (w).
*/


// NFFT output. Writes result from shared memory to global memory.

__device__ void 
NFFT_output_parallel( const unsigned int globalThreadId, const unsigned int domain_size_coils, const unsigned int number_of_samples, const unsigned int double_warp_size_power, const unsigned int sharedMemFirstSampleIdx, cuFloatComplex *out_data )
{
	// Output samples' complex values to global memory. Each thread writes 'domain_size_coils' complex floats.
	// All shared memory floats corresponding to domain 'threadIdx' is located in bank threadIdx%warp_size, i.e. no bank conflicts.

	for( unsigned int i=0; i<domain_size_coils; i++ ){
		cuFloatComplex sample_value = make_cuFloatComplex( shared_mem[sharedMemFirstSampleIdx+(i<<double_warp_size_power)], shared_mem[sharedMemFirstSampleIdx+(i<<double_warp_size_power)+my_warpSize] );
		out_data[i*number_of_samples+globalThreadId] = sample_value;
	}
}

__device__ void
__resolve_wrap( uint2 *grid_position, uint2 matrix_size_os )
{
	if( (int)(grid_position->x) < 0 )
		grid_position->x += matrix_size_os.x;

	if( (int)(grid_position->x) > matrix_size_os.x-1 )
		grid_position->x -= matrix_size_os.x;

	if( (int)(grid_position->y) < 0 )
		grid_position->y += matrix_size_os.y;

	if( (int)(grid_position->y) > matrix_size_os.y-1 )
		grid_position->y -= matrix_size_os.y;
}

__device__ void
__resolve_wrap( uint3 *grid_position, uint3 matrix_size_os )
{
	if( (int)(grid_position->x) < 0 )
		grid_position->x += matrix_size_os.x;

	if( (int)(grid_position->x) > matrix_size_os.x-1 )
		grid_position->x -= matrix_size_os.x;

	if( (int)(grid_position->y) < 0 )
		grid_position->y += matrix_size_os.y;

	if( (int)(grid_position->y) > matrix_size_os.y-1 )
		grid_position->y -= matrix_size_os.y;

	if( (int)(grid_position->z) < 0 )
		grid_position->z += matrix_size_os.z;

	if( (int)(grid_position->z) > matrix_size_os.z-1 )
		grid_position->z -= matrix_size_os.z;
}

__device__ void
__resolve_wrap( uint4 *grid_position, uint4 matrix_size_os )
{
	if( (int)(grid_position->x) < 0 )
		grid_position->x += matrix_size_os.x;

	if( (int)(grid_position->x) > matrix_size_os.x-1 )
		grid_position->x -= matrix_size_os.x;

	if( (int)(grid_position->y) < 0 )
		grid_position->y += matrix_size_os.y;

	if( (int)(grid_position->y) > matrix_size_os.y-1 )
		grid_position->y -= matrix_size_os.y;

	if( (int)(grid_position->z) < 0 )
		grid_position->z += matrix_size_os.z;

	if( (int)(grid_position->z) > matrix_size_os.z-1 )
		grid_position->z -= matrix_size_os.z;

	if( (int)(grid_position->w) < 0 )
		grid_position->w += matrix_size_os.w;

	if( (int)(grid_position->w) > matrix_size_os.w-1 )
		grid_position->w -= matrix_size_os.w;
}

// Parallel case, i.e. matrix_size_samples==1 && matrix_size_coils>1

template <class UINTd, class FLOATd, char TYPE, bool CORE> inline __device__ void
NFFT_parallel( const unsigned int globalThreadId, const unsigned int sharedMemFirstSampleIdx, const float alpha, const float W, const float beta, const unsigned int domain_size_coils, const UINTd matrix_size, const UINTd matrix_size_os, const unsigned int double_warp_size_power, const float one_over_num_encoding_steps, const unsigned int number_of_strips, const uint2 *stripsMap, const UINTd *stripsDir, const UINTd *stripOrigins, const unsigned int *stripLengths, bool lastDimFixed )
{
	const float one_over_W = 1.0f/W;
	const UINTd bias = matrix_size>>1;
	const UINTd bias_os = matrix_size_os>>1;
	const FLOATd matrix_size_os_f = uintd_to_floatd(matrix_size_os);

	// Sample position to convolve onto
	FLOATd sample_position;
	switch(TYPE){
		case 0:
			sample_position = NFFT_get_sample_position<FLOATd, CORE>( globalThreadId );
			break;
		case 1:
		case 2:
			sample_position = cuda_compute_radial_sample_position<UINTd, FLOATd, TYPE>( globalThreadId, alpha, bias, bias_os );
			break;
		default:
			return;
	}

	// Lookup the coordinate for the first sample strip for this domain (.x) and number of strips (.y)
	const uint2 strip_Index_Count = stripsMap[globalThreadId];

#ifdef BUFFER_SAFETY_CHECK
	if( (strip_Index_Count.y>0) && (((strip_Index_Count.x>=number_of_strips) || (strip_Index_Count.y>(UINT_MAX-strip_Index_Count.x)) || (strip_Index_Count.x+strip_Index_Count.y-1>=number_of_strips ))))
		return;
#endif

	// Some symbolic names
	const unsigned int firstStripIdx = strip_Index_Count.x;
	const unsigned int local_number_of_strips = strip_Index_Count.y;

	// Get direction for this thread's strips
	const UINTd strip_dir = stripsDir[globalThreadId];

#ifdef BUFFER_SAFETY_CHECK
	if( sum(strip_dir)!=1 )
		return;
#endif

	// Convolve each grid strip onto the domain samples (shared memory)
	for( unsigned int i=0; i<local_number_of_strips; i++ ){

		// Current strip index
		unsigned int current_strip_idx = firstStripIdx+i;
		
		// Get the current grid cell strip origin
		UINTd strip_origin = stripOrigins[current_strip_idx];

		// Get the current grid cell strip length
		unsigned int strip_length = stripLengths[current_strip_idx];

#ifdef BUFFER_SAFETY_CHECK
// What do about wrapping?
//		if( weak_greater_equal(strip_origin + (strip_length-1)*strip_dir, matrix_size_os ) )
//			return;
#endif

		// Iterate through all grid cells in current strip.
		for( unsigned int j=0; j<strip_length; j++ ){

			// Find grid position from the strip origin and direction.
			UINTd grid_position = strip_origin + j*strip_dir;

			// Grid position as FLOATd. (as _int_ to support wrapping)
			FLOATd float_grid_position = intd_to_floatd(uint_to_int(grid_position));

			// Resolve wrapping of grid position
			__resolve_wrap( &grid_position, matrix_size_os );

			// Calculate the distance between current sample and the grid cell
			FLOATd delta = fabsf(sample_position-float_grid_position);

			// If cell too distant from sample then move on to the next cell
			if( weak_greater(delta, W/2.0f ))
				continue;

			// Compute convolution weight.
			float weight = KaiserBessel( delta, matrix_size_os_f, one_over_W, beta, lastDimFixed );

			// Safety measure. We have occationally observed a NaN from the KaiserBessel computation
			if( !isfinite(weight) )
				continue;

			for( unsigned int image=0; image<domain_size_coils; image++ ){

				// Read the grid cell value from global memory
				cuFloatComplex grid_value = tex1Dfetch( TEX_GRID_VALUES_NFFT(CORE), image*prod(matrix_size_os)+co_to_idx( grid_position, matrix_size_os ) );

				// Add 'weight*sample_val' to the grid cell in shared memory (no bank conflicts)
				shared_mem[sharedMemFirstSampleIdx+(image<<double_warp_size_power)] += (weight*grid_value.x);
				shared_mem[sharedMemFirstSampleIdx+(image<<double_warp_size_power)+my_warpSize] += (weight*grid_value.y);
			}
		}
	}
}

template <char TYPE, bool CORE> inline __device__ void
NFFT_parallel_generic( const unsigned int globalThreadId, const unsigned int sharedMemFirstSampleIdx, const float alpha, const float W, const float beta, const unsigned int domain_size_coils, const uint2 matrix_size, const uint2 matrix_size_os, const uint2 matrix_size_wrap, const unsigned int double_warp_size_power, float2 *traj_positions, bool lastDimFixed )
{
	// Useful constants
	const float one_over_W = 1.0f/W;
	const float half_W = 0.5f*W;
	const float2 matrix_size_os_f = uintd_to_floatd(matrix_size_os);

	// Sample position to convolve onto
	// Computed in NFFT_H preprocessing, which includes a wrap zone. Remove this zone again...
	float2 sample_position = traj_positions[globalThreadId]-uintd_to_floatd(matrix_size_wrap>>1);

	int lower_limit_x = (int)ceil(sample_position.x-half_W);
	int lower_limit_y = (int)ceil(sample_position.y-half_W);
	int upper_limit_x = (int)floor(sample_position.x+half_W);
	int upper_limit_y = (int)floor(sample_position.y+half_W);

	// Iterate through all grid cells in current strip.
	for( int y=lower_limit_y; y<=upper_limit_y; y++ ){
		for( int x=lower_limit_x; x<=upper_limit_x; x++ ){

			// Find grid position from the strip origin and direction.
			uint2 grid_position = make_uint2( (unsigned int)x, (unsigned int)y );

			// Grid position as FLOATd. (as _int_ to support wrapping)
			float2 float_grid_position = intd_to_floatd(uint_to_int(grid_position));

			// Resolve wrapping of grid position
			__resolve_wrap( &grid_position, matrix_size_os );

			// Calculate the distance between current sample and the grid cell
			float2 delta = fabsf(sample_position-float_grid_position);

			// If cell too distant from sample then move on to the next cell
			if( weak_greater(delta, half_W ))
				continue;

			// Compute convolution weight.
			float weight = KaiserBessel( delta, matrix_size_os_f, one_over_W, beta, lastDimFixed );

			// Safety measure. We have occationally observed a NaN from the KaiserBessel computation
			if( !isfinite(weight) )
				continue;

			for( unsigned int image=0; image<domain_size_coils; image++ ){

				// Read the grid cell value from global memory
				cuFloatComplex grid_value = tex1Dfetch( TEX_GRID_VALUES_NFFT(CORE), image*prod(matrix_size_os)+co_to_idx( grid_position, matrix_size_os ) );

				// Add 'weight*sample_val' to the grid cell in shared memory (no bank conflicts)
				shared_mem[sharedMemFirstSampleIdx+(image<<double_warp_size_power)] += (weight*grid_value.x);
				shared_mem[sharedMemFirstSampleIdx+(image<<double_warp_size_power)+my_warpSize] += (weight*grid_value.y);
			}
		}
	}
}

template <char TYPE, bool CORE> inline __device__ void
NFFT_parallel_generic( const unsigned int globalThreadId, const unsigned int sharedMemFirstSampleIdx, const float alpha, const float W, const float beta, const unsigned int domain_size_coils, const uint3 matrix_size, const uint3 matrix_size_os, const uint3 matrix_size_wrap, const unsigned int double_warp_size_power, float3 *traj_positions, bool lastDimFixed )
{
}

template <char TYPE, bool CORE> inline __device__ void
NFFT_parallel_generic( const unsigned int globalThreadId, const unsigned int sharedMemFirstSampleIdx, const float alpha, const float W, const float beta, const unsigned int domain_size_coils, const uint4 matrix_size, const uint4 matrix_size_os, const uint4 matrix_size_wrap, const unsigned int double_warp_size_power, float4 *traj_positions, bool lastDimFixed )
{
}

/*
	These are the 'main' kernels for the NFFT
*/

// Generic kernel. 
// Special case where domain_size_samples is 1.

template <class UINTd, class FLOATd, char TYPE, bool CORE> __global__ void
NFFT_convolve_kernel_parallel_generic( float alpha, float W, float beta, unsigned int domain_count_samples, unsigned int domain_size_coils, UINTd matrix_size, UINTd matrix_size_os, UINTd matrix_size_wrap, unsigned int double_warp_size_power, FLOATd *traj_positions, cuFloatComplex *out_data, bool lastDimFixed )
{  
	// Global thread number	
	const unsigned int globalThreadId = (blockIdx.x*blockDim.x+threadIdx.x);

	// We might have a couple of warps left over due to our selection of grid/block dimensions.
	// This check will only work when there is NO __SYNCTHREADS() in the code!
	if( globalThreadId >= domain_count_samples )
		return;

	const unsigned int num_floats = domain_size_coils<<1;

	// All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size (i.e. no bank conflicts in warp).
	const unsigned int scatterSharedMemStart = (threadIdx.x/my_warpSize)*my_warpSize;
	const unsigned int scatterSharedMemStartOffset = threadIdx.x&(my_warpSize-1); //threadIdx.x%my_warpSize;
	const unsigned int sharedMemFirstSampleIdx = scatterSharedMemStart*num_floats + scatterSharedMemStartOffset;

	// Initialize shared memory
	for( unsigned int i=0; i<num_floats; i++ )
		shared_mem[sharedMemFirstSampleIdx+my_warpSize*i] = 0.0f;

	// Compute NFFT using arbitrary distributions (lookup trajectory positions in texture). The result goes to 'sample_values' in shared memory.
	NFFT_parallel_generic<TYPE, CORE>
		( globalThreadId, sharedMemFirstSampleIdx, alpha, W, beta, domain_size_coils, matrix_size, matrix_size_os, matrix_size_wrap, double_warp_size_power, traj_positions, lastDimFixed );

	// Output k-space image to global memory
	NFFT_output_parallel( globalThreadId, domain_size_coils, domain_count_samples, double_warp_size_power, sharedMemFirstSampleIdx, out_data );
}

// Parallel kernel. 
// Special case where domain_size_samples is 1.

template <class UINTd, class FLOATd, char TYPE, bool CORE> __global__ void
NFFT_convolve_kernel_parallel( float alpha, float W, float beta, unsigned int domain_count_samples, unsigned int domain_size_coils, UINTd matrix_size, UINTd matrix_size_os, unsigned int double_warp_size_power, float one_over_num_encodings, unsigned int number_of_strips, uint2 *stripsMap, UINTd *stripsDir, UINTd *stripOrigins, unsigned int *stripLengths, cuFloatComplex *out_data, bool lastDimFixed )
{  
	// Global thread number	
	const unsigned int globalThreadId = (blockIdx.x*blockDim.x+threadIdx.x);

	// We might have a couple of warps left over due to our selection of grid/block dimensions.
	// This check will only work when there is NO __SYNCTHREADS() in the code!
	if( globalThreadId >= domain_count_samples )
		return;

	const unsigned int num_floats = domain_size_coils<<1;

	// All shared memory floats corresponding to domain 'threadIdx.x' is located in bank threadIdx.x%warp_size (i.e. no bank conflicts in warp).
	const unsigned int scatterSharedMemStart = (threadIdx.x/my_warpSize)*my_warpSize;
	const unsigned int scatterSharedMemStartOffset = threadIdx.x&(my_warpSize-1); //threadIdx.x%my_warpSize;
	const unsigned int sharedMemFirstSampleIdx = scatterSharedMemStart*num_floats + scatterSharedMemStartOffset;

	// Initialize shared memory
	for( unsigned int i=0; i<num_floats; i++ )
		shared_mem[sharedMemFirstSampleIdx+my_warpSize*i] = 0.0f;

	// Compute NFFT using arbitrary distributions (lookup trajectory positions in texture). The result goes to 'sample_values' in shared memory.
	NFFT_parallel<UINTd, FLOATd, TYPE, CORE>
		( globalThreadId, sharedMemFirstSampleIdx, alpha, W, beta, domain_size_coils, matrix_size, matrix_size_os, double_warp_size_power, one_over_num_encodings, number_of_strips, stripsMap, stripsDir, stripOrigins, stripLengths, lastDimFixed );

	// Output k-space image to global memory
	NFFT_output_parallel( globalThreadId, domain_size_coils, domain_count_samples, double_warp_size_power, sharedMemFirstSampleIdx, out_data );
}
