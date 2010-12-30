#include "preprocess_radial.hcu"
#include "preprocess.hcu"
#include "preprocess_private.hcu"
#include "FFT.hcu"
#include "NFFT.hcu"
#include "NFFT_private.hcu"
#include "NSense_private.hcu"
#include "image_utilities.hcu"
#include "uint_util.hcu"
#include "uint_util_device.hcu"
#include "int_util_device.hcu"
#include "float_util.hcu"
#include "float_util_device.hcu"
#include "float_util_host.hcu"
#include "radixSort/radixsort.cuh"

#include <cublas.h>
#include <cudpp/cudpp.h>
#include <math_constants.h>
#include <vector_functions.h>

#include <math.h>
#include <assert.h>
#include <stdio.h>

const unsigned int preprocess_radial_NFFT_H_block_dim_x = 256;
CUDPPHandle scanplan;

/*
	Variables specific for the "online radial" types
*/

__constant__ float __num_samples_per_projection;
__constant__ float __num_samples_per_frame;
__constant__ float __num_projections_per_frame;
__constant__ float __one_over_num_samples_per_projection;
__constant__ float __one_over_num_samples_per_frame;
__constant__ float __one_over_num_projections_per_frame;
__constant__ float __angular_offset;
__constant__ float __gc_factor;
__constant__ float __frames_per_rotation_cycle;
__constant__ float __rotation_gap;
__constant__ float __interframe_rotation;
__constant__ float __one_over_radial_oversampling;


template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ void _preproc_NFFT_set_constants( PLAN<UINTd, FLOATd, TYPE> *plan )
{
	// Utility to set constant memory space

	float _one_over_radial_oversampling = (float)plan->matrix_size.x/(float)plan->samples_per_projection;
	cudaMemcpyToSymbol( __one_over_radial_oversampling, &_one_over_radial_oversampling, sizeof(float) );

	float _num_samples_per_projection = (float)plan->samples_per_projection;
	cudaMemcpyToSymbol( __num_samples_per_projection, &_num_samples_per_projection, sizeof(float) );

	float _num_samples_per_frame = (float)(plan->samples_per_projection*plan->projections_per_frame);
	cudaMemcpyToSymbol( __num_samples_per_frame, &_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_frame = 1.0f/(float)(plan->samples_per_projection*plan->projections_per_frame);
	cudaMemcpyToSymbol( __one_over_num_samples_per_frame, &_one_over_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_projection = 1.0f/(float)plan->samples_per_projection;
	cudaMemcpyToSymbol( __one_over_num_samples_per_projection, &_one_over_num_samples_per_projection, sizeof(float) );

	float _one_over_num_projections_per_frame = 1.0f/(float)plan->projections_per_frame;
	cudaMemcpyToSymbol( __one_over_num_projections_per_frame, &_one_over_num_projections_per_frame, sizeof(float) );

	float _num_projections_per_frame = (float)(plan->projections_per_frame);
	cudaMemcpyToSymbol( __num_projections_per_frame, &_num_projections_per_frame, sizeof(float) );

	switch(TYPE){
	
		case 1:
			{
				cudaMemcpyToSymbol( __gc_factor, &plan->gc_factor, sizeof(float) );

				float _angular_offset = (float)plan->angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		case 2:
			{
				cudaMemcpyToSymbol( __interframe_rotation, &plan->interframe_rotation, sizeof(float) );

				float _frames_per_rotation_cycle = (float)plan->frames_per_rotation;
				cudaMemcpyToSymbol( __frames_per_rotation_cycle, &_frames_per_rotation_cycle, sizeof(float) );

				float _rotation_gap = (float)((plan->frames_per_rotation-1)/2);
				cudaMemcpyToSymbol( __rotation_gap, &_rotation_gap, sizeof(float) );

				float _angular_offset = (float)plan->angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		default:
			printf("\nInternal ERROR!: 'NFFT_set_constants': switch. Quitting.\n" ); fflush(stdout); 
			exit(1);
			break;
	}	

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nInternal ERROR!: 'NFFT_set_constants': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}
}

template< char TYPE > 
__host__ void _preproc_NFFT_set_constants( unsigned int matrix_size, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset = 0, unsigned int frames_per_rotation = 1, float gc_factor = 1.0f )
{
	// Utility to set constant memory space

	float _one_over_radial_oversampling = (float)matrix_size/(float)samples_per_projection;
	cudaMemcpyToSymbol( __one_over_radial_oversampling, &_one_over_radial_oversampling, sizeof(float) );

	float _num_samples_per_projection = (float)samples_per_projection;
	cudaMemcpyToSymbol( __num_samples_per_projection, &_num_samples_per_projection, sizeof(float) );

	float _num_samples_per_frame = (float)(samples_per_projection*projections_per_frame);
	cudaMemcpyToSymbol( __num_samples_per_frame, &_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_frame = 1.0f/(float)_num_samples_per_frame;
	cudaMemcpyToSymbol( __one_over_num_samples_per_frame, &_one_over_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_projection = 1.0f/(float)samples_per_projection;
	cudaMemcpyToSymbol( __one_over_num_samples_per_projection, &_one_over_num_samples_per_projection, sizeof(float) );

	float _one_over_num_projections_per_frame = 1.0f/(float)projections_per_frame;
	cudaMemcpyToSymbol( __one_over_num_projections_per_frame, &_one_over_num_projections_per_frame, sizeof(float) );

	float _num_projections_per_frame = (float)(projections_per_frame);
	cudaMemcpyToSymbol( __num_projections_per_frame, &_num_projections_per_frame, sizeof(float) );

	switch(TYPE){
	
		case 1:
			{
				cudaMemcpyToSymbol( __gc_factor, &gc_factor, sizeof(float) );

				float _angular_offset = (float)angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		case 2:
			{
				float _interframe_rotation = CUDART_PI_F/(float)(frames_per_rotation*projections_per_frame);
				cudaMemcpyToSymbol( __interframe_rotation, &_interframe_rotation, sizeof(float) );

				float _frames_per_rotation_cycle = (float)frames_per_rotation;
				cudaMemcpyToSymbol( __frames_per_rotation_cycle, &_frames_per_rotation_cycle, sizeof(float) );

				float _rotation_gap = (float)((frames_per_rotation-1)/2);
				cudaMemcpyToSymbol( __rotation_gap, &_rotation_gap, sizeof(float) );

				float _angular_offset = (float)angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		default:
			printf("\nInternal ERROR!: 'NFFT_set_constants': switch. Quitting.\n" ); fflush(stdout); 
			exit(1);
			break;
	}	

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nInternal ERROR!: 'NFFT_set_constants': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}
}

template< class UINTd, class FLOATd, char TYPE > mr_recon::NFFT_plan<UINTd, FLOATd, TYPE>*
preprocess_radial_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	if( sizeof(UINTd)==sizeof(uint2) && projections_per_frame != num_projections ){
		printf("\nError: 'preprocess_radial_NFFT': For 2D reconstructions the projections must reflect a single frame.\n");
		return 0x0;
	}

	mr_recon::NFFT_plan<UINTd, FLOATd, TYPE> *tmp = new mr_recon::NFFT_plan<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_samples, domain_size_coils, W );

	tmp->preprocess_radial( num_projections, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}

template< class UINTd, class FLOATd, char TYPE > mr_recon::NFFT_H_plan<UINTd, FLOATd, TYPE>*
preprocess_radial_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_coils, float W, unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{	
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	if( sizeof(UINTd)==sizeof(uint2) && projections_per_frame != num_projections ){
		printf("\nError: 'preprocess_radial_NFFT': For 2D reconstructions the projections must reflect a single frame.\n");
		return 0x0;
	}

	mr_recon::NFFT_H_plan<UINTd, FLOATd, TYPE> *tmp = new mr_recon::NFFT_H_plan<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_coils, W );

	tmp->preprocess_radial( num_projections, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}

template< class UINTd, class FLOATd, char TYPE > mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE>*
preprocess_radial_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	if( sizeof(UINTd)==sizeof(uint2) && projections_per_frame != num_projections ){
		printf("\nError: 'preprocess_radial_NFFT': For 2D reconstructions the projections must reflect a single frame.\n");
		return 0x0;
	}

	mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *tmp = new mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_samples, domain_size_coils, W );

	tmp->preprocess_radial( num_projections, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}


// Utilities for 'preprocess_radial_trajectories_NFFT_kernel'
// Notice: Hard coded to the current assumptions about fixed dimensions!

__inline__ __device__ uint2 __make_bias( uint2 size )
{
	return make_uint2( size.x>>1, size.y>>1 );
}

__inline__ __device__ uint3 __make_bias( uint3 size )
{
	return make_uint3( size.x>>1, size.y>>1, 0 );
}

__inline__ __device__ void __make_x_vec( uint2 *x )
{
	*x = make_uint2( 1,0 );
}

__inline__ __device__ void __make_x_vec( uint3 *x )
{
	*x = make_uint3( 1,0,0 );
}

__inline__ __device__ void __make_y_vec( uint2 *y )
{
	*y = make_uint2( 0,1 );
}

__inline__ __device__ void __make_y_vec( uint3 *y )
{
	*y = make_uint3( 0,1,0 );
}

__inline__ __device__ void __set_origin( uint2 *ptr, unsigned int x, unsigned int y, unsigned int z )
{
	*ptr = make_uint2( x, y );
}

__inline__ __device__ void __set_origin( uint3 *ptr, unsigned int x, unsigned int y, unsigned int z )
{
	*ptr = make_uint3( x, y, z );
}

template <class UINTd> inline __device__ 
bool isXdir( UINTd stripDir )
{
	return (stripDir.x == 1);
}

template <class UINTd> inline __device__ 
bool isYdir( UINTd stripDir )
{
	return (stripDir.y == 1);
}

template< class UINTd, class FLOATd, char TYPE > __global__ void
preprocess_radial_trajectories_NFFT_kernel( unsigned int max_strips_per_domain, unsigned int domain_count, unsigned int domain_size, UINTd matrix_size, UINTd matrix_size_os, float alpha, float kernel_width, uint2 *stripsMapPtr, UINTd *stripsDirPtr, UINTd *stripOriginsPtr, unsigned int *stripLengthsPtr )
{
	// Global thread number	
	const unsigned int globalThreadId = (blockIdx.x*blockDim.x+threadIdx.x);

	// We might have a couple of warps left over due to our selection of grid/block dimensions.
	// This check will only work when there is NO __SYNCTHREADS() in the code!
	if( globalThreadId >= domain_count )
		return;

	/*
		The algorithm is as follows
		---------------------------

		1.	Create bounding box for each "sample domain" based on the first and last sample extended by the haæf convolution kernel size
			All grid cells that influence the samples in the current domain is hereby contained in the bounding box.

		2.	Set the strip direction to the longest axis of the bounding box
			I.e. the number of strips will relate to the narrow dimension

		3. Iterate over rows/columns in the narrow dimension

			3.1.	Find the 'k-space projection vector p_r' between the first and last samples (x1,y1) and (x2,y2).
					Set v'=[y2-y1, -(x2-x1)], v = v'/|v| (perpendicular to 'p_r'), and q=(x0,y0) the start/end point (to be determined) on the strip. Either x0 or y0 is given.
					The start/end point is given by solving "a(x0+W/2*v.x)+b=y0+W/2*v.y" - i.e. which point 'q' is moved onto to the projection when moved W/2 along 'v'. 'a' and 'b' correspond to the line equation ax=b for the projection.

			3.2		Crop the two points to the bounsing box.

			3.3		Construct strip based on the two (x,y) pairs found in 3.1. Write to global memory.
		
		4.	Write "per-thread" datastructure to global memory.

	*/


	// Find first and last sample in projection section
	const unsigned int sampleIdx = globalThreadId*domain_size;
	const unsigned int num_samples_per_projection = float2uint(__num_samples_per_projection);
	const unsigned int num_projections_per_frame = float2uint(__num_projections_per_frame);

	unsigned int current_projection = sampleIdx/num_samples_per_projection;
	const unsigned int current_frame = current_projection/num_projections_per_frame;
	const UINTd bias = __make_bias(matrix_size);
	const float half_kernel_width = kernel_width/2.0f;
	float cos_angle, sin_angle;

	switch(TYPE){
	case 1:
		{
			const float angle_step = CUDART_PI_F/((sqrtf(5.0f) + 1.0f)/2.0f);
			__sincosf( (uint2float(current_projection)+__angular_offset)*angle_step*__gc_factor, &sin_angle, &cos_angle );
		}
		break;
	case 2:
		{
			const float cur_proj = uint2float(current_projection);
			const float frame = floorf((cur_proj+0.5f)*__one_over_num_projections_per_frame);
			const float cur_proj_local = cur_proj - frame*__num_projections_per_frame;
			const float rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
			__sincosf( cur_proj_local*CUDART_PI_F*__one_over_num_projections_per_frame+rotation*__interframe_rotation, &sin_angle, &cos_angle );
		}
		break;
	default:
		break;
	}

	const float2 samplePos_1 = make_float2( (((sampleIdx%num_samples_per_projection)*__one_over_radial_oversampling)-bias.x)*cos_angle*alpha, (((sampleIdx%num_samples_per_projection)*__one_over_radial_oversampling)-bias.y)*sin_angle*alpha );
	const float2 samplePos_n = make_float2( ((((sampleIdx+domain_size-1)%num_samples_per_projection)*__one_over_radial_oversampling)-bias.x)*cos_angle*alpha, ((((sampleIdx+domain_size-1)%num_samples_per_projection)*__one_over_radial_oversampling)-bias.y)*sin_angle*alpha );

	// 1. Define bounding box corners (origin is center of grid)
	float2 _bbox_1 = samplePos_1 + make_float2((samplePos_1.x<=samplePos_n.x) ? -half_kernel_width : half_kernel_width, (samplePos_1.y<=samplePos_n.y) ? -half_kernel_width : half_kernel_width);
	float2 _bbox_2 = samplePos_n + make_float2((samplePos_1.x<=samplePos_n.x) ? half_kernel_width : -half_kernel_width, (samplePos_1.y<=samplePos_n.y) ? half_kernel_width : -half_kernel_width);

	// confine _bbox to grid
	if(_bbox_1.x<_bbox_2.x){
		_bbox_1.x = ceil(_bbox_1.x);
		_bbox_2.x = floor(_bbox_2.x);
	}
	else{
		_bbox_1.x = floor(_bbox_1.x);
		_bbox_2.x = ceil(_bbox_2.x);
	}
	if(_bbox_1.y<_bbox_2.y){
		_bbox_1.y = ceil(_bbox_1.y);
		_bbox_2.y = floor(_bbox_2.y);
	}
	else{
		_bbox_1.y = floor(_bbox_1.y);
		_bbox_2.y = ceil(_bbox_2.y);
	}

	// Bias _bbox to [0;matrix_size_os] (plus wrapping)
	const FLOATd matrix_size_os_f = uintd_to_floatd(matrix_size_os);
	float2 bias_os_f = uintd_to_floatd(uintd_to_uint2(matrix_size_os>>1));
	_bbox_1 += bias_os_f;
	_bbox_2 += bias_os_f;

	// 2. Define strip direction to "match" the slope of the projection

	UINTd stripDir;

	// cos(PI/4)== sqrt(0.5)
	if( fabsf(cos_angle)>CUDART_SQRT_HALF_F )
		__make_x_vec( &stripDir );
	else
		__make_y_vec( &stripDir );

	// Write stripsDir to global memory
	stripsDirPtr[globalThreadId] = stripDir;

	// Make sure the bbox is not "inverted" in the strip direction
	const float2 _diff = _bbox_2-_bbox_1;
	const float2 _adiff = make_float2(fabsf(_diff.x), fabsf(_diff.y));

	if( isXdir(stripDir) ){
		if(_diff.x<0){
			float2 tmp = _bbox_1;
			_bbox_1 = _bbox_2;
			_bbox_2 = tmp;
		}
	}
	else{
		if(_diff.y<0){
			float2 tmp = _bbox_1;
			_bbox_1 = _bbox_2;
			_bbox_2 = tmp;
		}
	}

	// Cast _bbox to int2
	int2 bbox_1 = floatd_to_intd(_bbox_1);
	int2 bbox_2 = floatd_to_intd(_bbox_2);

	// Define v (perpendicular to projection)
	const float2 v = make_float2( sin_angle, -cos_angle );

	// 3. Iterate over the narrow dimension

	const unsigned int numIterations = (isXdir(stripDir)) ? float2uint(_adiff.y)+1 : float2uint(_adiff.x)+1;
	unsigned int num_strips = 0;

	for( int i=0; i<numIterations; i++ ){

		// 3.1. Find first/last grid cell for the strip

		if( isXdir(stripDir) ){
	
			// NOTE: uints are cast to int in convolution kernel

			const int y0 = min(bbox_1.y,bbox_2.y)+i; // min neccessary - ordered according to 'x'
			__set_origin( &stripOriginsPtr[globalThreadId*max_strips_per_domain+num_strips], (unsigned int)bbox_1.x, (unsigned int)y0, current_frame );
			stripLengthsPtr[globalThreadId*max_strips_per_domain+num_strips] = bbox_2.x-bbox_1.x+1;
			num_strips++;
		}

		else{

			// NOTE: uints are cast to int in convolution kernel

			const int x0 = min(bbox_1.x,bbox_2.x)+i; // min neccessary - ordered according to 'y'
			__set_origin( &stripOriginsPtr[globalThreadId*max_strips_per_domain+num_strips], (unsigned int)x0, (unsigned int)bbox_1.y, current_frame );
			stripLengthsPtr[globalThreadId*max_strips_per_domain+num_strips] = bbox_2.y-bbox_1.y+1;
			num_strips++;
		}
	}

	// Write num_strips to global memory
	stripsMapPtr[globalThreadId] = make_uint2( globalThreadId*max_strips_per_domain, num_strips );
}


template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_plan< UINTd, FLOATd, TYPE >::preprocess_radial( unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	successfully_preprocessed = false;

	number_of_samples = num_projections*samples_per_projection;
	domain_count_samples = number_of_samples/domain_size_samples;

	// Some sanity checks

	if( sizeof(UINTd)==sizeof(uint2) && projections_per_frame != num_projections ){
		printf("\nError: 'preprocess_radial': For 2D reconstructions the projections must reflect a single frame.\n");
		return 0x0;
	}

	if( sizeof(UINTd)<sizeof(uint3) && sum(fixed_dims) ){
		printf("\nERROR: NFFT_plan< UINTd, FLOATd, TYPE >::preprocess_radial : There can be no fixed dimensions in the 2D case.");
		return false;				
	}

	if( sizeof(UINTd)==sizeof(uint3) && !get_last_dim(fixed_dims) ){
		printf("\nERROR: NFFT_plan< UINTd, FLOATd, TYPE >::preprocess_radial : For now, the last dimension must be fixed in the 3D case.");
		return false;				
	}

	if( sum(fixed_dims)>1 || (sum(fixed_dims)==1 && !get_last_dim(fixed_dims)) ){
		printf("\nERROR: NFFT_plan< UINTd, FLOATd, TYPE >::preprocess_radial : Only one fixed dimension allowed. It must be the last dimension.");
		return false;				
	}

	if( number_of_samples%domain_size_samples ) {
		printf("\nERROR: We currently require the number of samples to be a multiple of the domain length.\n");
		return false;				
	}

	if( (number_of_samples/domain_size_samples)%2 ){
		printf("\nERROR: We currently require the number of samples to be an EVEN multiple of the domain length.\n");
		return false;				
	}

	// Update plan
	this->samples_per_projection = samples_per_projection;
	this->projections_per_frame = projections_per_frame;
	this->angular_offset = angular_offset;
	this->frames_per_rotation = frames_per_rotation;
	this->interframe_rotation = CUDART_PI_F/(float)(frames_per_rotation*projections_per_frame);
	this->gc_factor = gc_factor;
	this->total_projections_f = (float)num_projections;

	/*
		The maximum number of strips occurs for projections in a 45 degree angle.
		The projection causes a maximum of 1+ceil(sqrt((l^2)/2)) strips where l is the length of the projection section.
		The convolution kernel width account for a maximum of 2*floor(W/2) additional strips. 
	*/

	float l = (domain_size_samples-1.0f)*(matrix_size_os.x/(float)samples_per_projection);
	this->max_strips_per_domain_NFFT = 1 + (unsigned int)(ceil(sqrt((l*l)/2.0f))) + (unsigned int)floor(W);

	this->number_of_strips_NFFT = domain_count_samples*max_strips_per_domain_NFFT;

	assert((number_of_samples%domain_size_samples)==0);
	unsigned int number_of_threads = number_of_samples/domain_size_samples;

	if( stripsMapDevPtr_NFFT ){
		cudaFree( stripsMapDevPtr_NFFT );
		stripsMapDevPtr_NFFT = 0x0;
	}
	if( stripsDirDevPtr_NFFT ){
		cudaFree( stripsDirDevPtr_NFFT );
		stripsDirDevPtr_NFFT = 0x0;
	}
	if( stripOriginsDevPtr_NFFT ){
		cudaFree( stripOriginsDevPtr_NFFT );
		stripOriginsDevPtr_NFFT = 0x0;
	}
	if( stripLengthsDevPtr_NFFT ){
		cudaFree( stripLengthsDevPtr_NFFT );
		stripLengthsDevPtr_NFFT = 0x0;
	}

	cudaMalloc( (void**) &stripsMapDevPtr_NFFT, number_of_threads*sizeof(uint2) );
	cudaMalloc( (void**) &stripsDirDevPtr_NFFT, number_of_threads*sizeof(UINTd) );
	cudaMalloc( (void**) &stripOriginsDevPtr_NFFT, number_of_threads*max_strips_per_domain_NFFT*sizeof(UINTd));
	cudaMalloc( (void**) &stripLengthsDevPtr_NFFT, number_of_threads*max_strips_per_domain_NFFT*sizeof(unsigned int));

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error (malloc?) in 'preprocess_radial_trajectories_NFFT'!"); fflush(stdout);
		return false;
	}

	// Invoke kernel

	_preproc_NFFT_set_constants( this );

	unsigned int block_dim_x = 320;
	dim3 blockDim(block_dim_x,1,1);
	dim3 gridDim((unsigned int)ceil((double)number_of_threads/(double)block_dim_x), 1, 1 );

	preprocess_radial_trajectories_NFFT_kernel<UINTd,FLOATd,TYPE><<<gridDim, blockDim>>>( max_strips_per_domain_NFFT, number_of_threads, domain_size_samples, matrix_size, matrix_size_os, alpha, W, stripsMapDevPtr_NFFT, stripsDirDevPtr_NFFT, stripOriginsDevPtr_NFFT, stripLengthsDevPtr_NFFT );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'preprocess_radial_trajectories_NFFT_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	successfully_preprocessed = true;

	return true;
}

// Utilities for 'preprocess_radial_trajectories_NFFT_H_kernel'
// Notice: Hard coded to the current assumptions about fixed dimensions!

__inline__ __device__ void __set_bbbox( float2 *ptr, float x, float y, float z )
{
	*ptr = make_float2( x, y );
}

__inline__ __device__ void __set_bbbox( float3 *ptr, float x, float y, float z )
{
	*ptr = make_float3( x, y, z );
}

#define _BR(i) _br[(i)*preprocess_radial_NFFT_H_block_dim_x+threadIdx.x]
#define _R(i) _r[(i)*preprocess_radial_NFFT_H_block_dim_x+threadIdx.x]
#define R(i) r[(i)*preprocess_radial_NFFT_H_block_dim_x+threadIdx.x]

template< class UINTd, class FLOATd, char TYPE> __global__ void
preprocess_radial_trajectories_NFFT_H_kernel( UINTd matrix_size_os, UINTd domain_count_grid, UINTd domain_size_grid, UINTd fixed_dims, float kernel_width, unsigned int *stripOffsets, uint2 *stripsMapPtr, uint2 *stripsPtr )
{

	// Global thread number	
	const unsigned int globalThreadId = blockIdx.y*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x;

	// Number of domains
	const unsigned int number_of_domains = prod(domain_count_grid);

	// We might have a couple of warps left over due to our selection of grid/block dimensions.
	// This check will only work when there is NO __SYNCTHREADS() in the code!
	if( globalThreadId >= number_of_domains )
		return;

	// Find the global grid coordinate (upper left corner) of the thread's domain
	const UINTd _globalIdx = idx_to_co(globalThreadId, domain_count_grid);
	const unsigned int num_projections_per_frame = float2uint(__num_projections_per_frame);
	const unsigned int num_samples_per_projection = float2uint(__num_samples_per_projection);

	// Find domain position (make sure an enventual time dimension is unchanged)
	FLOATd globalIdx = uintd_to_floatd(_globalIdx*domain_size_grid);
	globalIdx.x -= uint2float((domain_count_grid.x*domain_size_grid.x)>>1);
	globalIdx.y -= uint2float((domain_count_grid.y*domain_size_grid.y)>>1);
	
	FLOATd bbox_ll, bbox_lr, bbox_ul, bbox_ur;

	__set_bbbox( &bbox_ll, globalIdx.x-kernel_width/2.0f, globalIdx.y-kernel_width/2.0f, get_last_dim(globalIdx) );
	__set_bbbox( &bbox_lr, globalIdx.x+domain_size_grid.x-1+kernel_width/2.0f, globalIdx.y-kernel_width/2.0f, get_last_dim(globalIdx) );
	__set_bbbox( &bbox_ul, globalIdx.x-kernel_width/2.0f, globalIdx.y+domain_size_grid.y-1+kernel_width/2.0f, get_last_dim(globalIdx) );
	__set_bbbox( &bbox_ur, globalIdx.x+domain_size_grid.x-1+kernel_width/2.0f, globalIdx.y+domain_size_grid.y-1+kernel_width/2.0f,	get_last_dim(globalIdx) );

//	printf("\nbbox_ul: %f %f %f", bbox_ul.x, bbox_ul.y, bbox_ul.z );
//	printf("\nbbox_lr: %f %f %f", bbox_lr.x, bbox_lr.y, bbox_lr.z );

	unsigned char scenario;
	float theta_min, theta_max;

	// Take special care if any of the bbox corners is exactly the origin (atan2f undefined then)
	if( (bbox_ul.x == 0.0f && bbox_ul.y == 0.0f) || (bbox_ur.x == 0.0f && bbox_ur.y == 0.0f) || (bbox_ll.x == 0.0f && bbox_ll.y == 0.0f) || (bbox_lr.x == 0.0f && bbox_lr.y == 0.0f) )
	{
		// Treat as if the origin is fully inside the bbox. All projections are hit.
		scenario = 1; 
	}
	else{

		// Find the angles corresponding to the bounding box corners
		float _theta1 = atan2f(bbox_ul.y,bbox_ul.x);
		float _theta2 = atan2f(bbox_ur.y,bbox_ur.x);
		float _theta3 = atan2f(bbox_ll.y,bbox_ll.x);
		float _theta4 = atan2f(bbox_lr.y,bbox_lr.x);

		//	printf("\ntheta: %f %f %f %f", _theta1, _theta2, _theta3, _theta4 );

		// Find minimum and maximum angles
		float _theta_min1, _theta_min2, _theta_max1, _theta_max2;

		_theta_min1 = min( _theta1, _theta2 );
		_theta_min2 = min( _theta3, _theta4 );
		_theta_max1 = max( _theta1, _theta2 );
		_theta_max2 = max( _theta3, _theta4 );

		theta_min = min( _theta_min1, _theta_min2 );
		theta_max = max( _theta_max1, _theta_max2 );

		if( ((theta_min*theta_max)>=0) || (theta_min>-CUDART_PIO2_F)){

			// Normal case, no discontinuity 
			scenario = 0;
		}
		else if(bbox_ur.x>=0){ 

			// bbox crosses both 'x' and 'y' axis, is thus "around" the origin
			scenario = 1;
		}
		else{

			// bbox crosses 'x' axis with all negative 'x' values. Thus 'theta' has a discontinuity.
			scenario = 2;

			if(_theta1<0)
				_theta1 += (2.0f*CUDART_PI_F);
			if(_theta2<0)
				_theta2 += (2.0f*CUDART_PI_F);
			if(_theta3<0)
				_theta3 += (2.0f*CUDART_PI_F);
			if(_theta4<0)
				_theta4 += (2.0f*CUDART_PI_F);

			_theta_min1 = min( _theta1, _theta2 );
			_theta_min2 = min( _theta3, _theta4 );
			_theta_max1 = max( _theta1, _theta2 );
			_theta_max2 = max( _theta3, _theta4 );

			theta_min = min( _theta_min1, _theta_min2 );
			theta_max = max( _theta_max1, _theta_max2 );
		}
	}

	unsigned int cur_projection = get_last_dim(fixed_dims)*get_last_dim(_globalIdx)*num_projections_per_frame;
	unsigned int num_strips = 0;

	// Test all projections for intersection with bounding box
	while( cur_projection < num_projections_per_frame*(get_last_dim(fixed_dims)*get_last_dim(_globalIdx)+1) ){

		bool found = false;
		float theta;

		// Find first/next projection inside bounding box
		for( ; cur_projection<num_projections_per_frame*(get_last_dim(fixed_dims)*get_last_dim(_globalIdx)+1); cur_projection++ ){

			// Calculate angle for current projection [-PI;PI[

			switch(TYPE){
			case 1:
				{
					const float angle_step = CUDART_PI_F/((sqrtf(5.0f) + 1.0f)/2.0f); // Golden angle step size
					theta = fmodf( (uint2float(cur_projection)+__angular_offset)*angle_step*__gc_factor, 2.0f*CUDART_PI_F );
				}
				break;

			case 2:
				{
					const float cur_proj = uint2float(cur_projection);
					const float frame = floorf((cur_proj+0.5f)*__one_over_num_projections_per_frame);
					const float cur_proj_local = cur_proj - frame*__num_projections_per_frame;
					const float rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
					theta = fmodf( cur_proj_local*CUDART_PI_F*__one_over_num_projections_per_frame+rotation*__interframe_rotation, 2.0f*CUDART_PI_F );
				}
				break;

			default:
				break;
			}

			if( theta>CUDART_PI_F )
				theta -= (2.0f*CUDART_PI_F);

			// Check if theta is inside bounding box (plus mirrored box)
			switch( scenario ){

				case 0:

					// Normal case, no discontinuity 
					if( (theta>theta_min && theta<theta_max) || (theta>(theta_min+CUDART_PI_F) && theta<(theta_max+CUDART_PI_F)) || (theta>(theta_min-CUDART_PI_F) && theta<(theta_max-CUDART_PI_F)) )
						found = true;
					break;

				case 1:

					// bbox crosses both 'x' and 'y' axis, is thus "around" the origin
					found = true;
					break;
				
				case 2:

					// bbox crosses 'x' axis with all negative 'x' values. Thus 'theta' has a discontinuity.
					if(theta<0)
						theta += (2.0f*CUDART_PI_F);

					// Recheck if theta is inside bounding box (plus mirrored box). Theta in [0;2*PI[
					if( (theta>theta_min && theta<theta_max) || (theta>(theta_min+CUDART_PI_F) && theta<(theta_max+CUDART_PI_F)) || (theta>(theta_min-CUDART_PI_F) && theta<(theta_max-CUDART_PI_F)) )
						found = true;
					break;
			
				default:
					break;
			}

			// When the first intersecting projection is found we break out the for-loop.
			if( found )
				break;
		}

		if( found ){ 
			
			// Projection enters bounding box

			/*
				Find radii for the four projection-line intersections of the bbox:
					- intersection with vertical line: r*cos(theta)=x
					- intersection with horizontal line: r*sin(theta)=y
			*/

			float _cos_angle, _sin_angle;
			__sincosf( theta, &_sin_angle, &_cos_angle );

			const float __limit = 0.000001f;

			__shared__ bool _br[4*preprocess_radial_NFFT_H_block_dim_x]; // Accessed using the _BR macro

			if( fabsf(_cos_angle)<__limit )
				_BR(0) = _BR(1) = false;
			else
				_BR(0) = _BR(1) = true;

			if( fabsf(_sin_angle)<__limit )
				_BR(2) = _BR(3) = false;
			else
				_BR(2) = _BR(3) = true;

			__shared__ float _r[4*preprocess_radial_NFFT_H_block_dim_x]; // Accessed using the _R macro

			if( _BR(0) )
				_R(0) = bbox_ll.x / _cos_angle;
			if( _BR(1) )
				_R(1) = bbox_ur.x / _cos_angle;
			if( _BR(2) )
				_R(2) = bbox_ll.y / _sin_angle;
			if( _BR(3) )
				_R(3) = bbox_ur.y / _sin_angle;

			// Determine which radii are within the allowed interval?
			if(_BR(0)) 
				_BR(0) = (_R(0)*_sin_angle>=bbox_ll.y) && (_R(0)*_sin_angle<=bbox_ul.y);
			if(_BR(1)) 
				_BR(1) = (_R(1)*_sin_angle>=bbox_lr.y) && (_R(1)*_sin_angle<=bbox_ur.y);
			if(_BR(2)) 
				_BR(2) = (_R(2)*_cos_angle>=bbox_ll.x) && (_R(2)*_cos_angle<=bbox_lr.x);
			if(_BR(3)) 
				_BR(3) = (_R(3)*_cos_angle>=bbox_ul.x) && (_R(3)*_cos_angle<=bbox_ur.x);

// Emulation mode only
//			assert( (_BR(0) + _BR(1) + _BR(2) + _BR(3)) == 2 );
//			if( (_BR(0) + _BR(1) + _BR(2) + _BR(3)) != 2 )
//				_BR(0)!=_BR(0); // ERROR STATE!!!
// end: emulation mode only

			// Radii of intersection points
			__shared__ float r[2*preprocess_radial_NFFT_H_block_dim_x];	// Accessed using the R macro

			unsigned int hits = _BR(0) + _BR(1) +_BR(2) + _BR(3);

			// Find r's and make sure r[0]<r[1]
			if( hits < 2 ){
				// Something is wrong if we get in here!!!
				cur_projection++;
				continue;
			}
			else{
				unsigned int found2 = 0;
				for( unsigned int i=0; i<4; i++ ){
					if( _BR(i) ){
						if( found2 == 0 ){
							R(0) = _R(i);
						}
						else if( found2 == 1 ){
							if( _R(i)<R(0) ){
								R(1) = R(0);
								R(0) = _R(i);
							}
							else{
								R(1) = _R(i);
							}
						}
						else{
							if( _R(i)<R(0) )
								R(0) = _R(i);
							else if( _R(i)>R(1) )
								R(1) = _R(i);
						}
						
						found2++;
					}
				}
			}

			// Convert the two radii in 'r' to sample ids		
			// left in sum: midpoint idx. right in sum: signed offset from midpoint.
			int sample_idx[2];
			sample_idx[0] = (num_samples_per_projection>>1)+float2int(ceil(R(0)/(uint2float(matrix_size_os.x)/__num_samples_per_projection)));
			sample_idx[1] = (num_samples_per_projection>>1)+float2int(floor(R(1)/(uint2float(matrix_size_os.x)/__num_samples_per_projection)));

			// If sample_idx[0]>sample_idx[1] the strip is empty (can happen due to the ceil/floor)
			unsigned int outside = 0 + (sample_idx[0]>sample_idx[1])*2;

			// Make sure sample_idx falls in the allowed range
			for( unsigned int i=0; i<2; i++ ){

				if( sample_idx[i]<0 ){
					sample_idx[i]=0;
					outside++;
				}

				if( sample_idx[i]>(num_samples_per_projection-1) ){
					sample_idx[i]=num_samples_per_projection-1;
					outside++;
				}
			}

			// If both samples are outside the valid range we discard this projection (we assume that no kernel will cover the whole image)
			if( outside < 2 ){

				// Convert sample ids to strip entry
				// We will use global indices for now for compatibility with "offline code".
				uint2 strip = make_uint2( cur_projection*num_samples_per_projection+sample_idx[0], sample_idx[1]-sample_idx[0]+1 );

				// Write strip to global memory (using coalesced write)
				stripsPtr[stripOffsets[globalThreadId]+num_strips] = strip;
				num_strips++;
			}
			cur_projection++;
		}
	}

	// Write num_strips to global memory
	stripsMapPtr[globalThreadId] = make_uint2( stripOffsets[globalThreadId], num_strips );
}


template< class UINTd, class FLOATd, char TYPE> __global__ void
preprocess_radial_allocate_kernel( UINTd matrix_size_os, UINTd domain_count_grid, UINTd domain_size_grid, UINTd fixed_dims, float kernel_width, unsigned int *stripCounts, uint2 *stripCountPairs )
{
	// Code duplicated from preprocess_radial kernel

	// Global thread number	
	const unsigned int globalThreadId = blockIdx.y*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x;

	// Number of domains
	const unsigned int number_of_domains = prod(domain_count_grid);

	// We might have a couple of warps left over due to our selection of grid/block dimensions.
	// This check will only work when there is NO __SYNCTHREADS() in the code!
	if( globalThreadId >= number_of_domains )
		return;

	// Find the global grid coordinate (upper left corner) of the thread's domain
	const UINTd _globalIdx = idx_to_co(globalThreadId, domain_count_grid);
	const unsigned int num_projections_per_frame = float2uint(__num_projections_per_frame);
	const unsigned int num_samples_per_projection = float2uint(__num_samples_per_projection);

	// Find domain position (make sure an enventual time dimension is unchanged)
	FLOATd globalIdx = uintd_to_floatd(_globalIdx*domain_size_grid);
	globalIdx.x -= uint2float((domain_count_grid.x*domain_size_grid.x)>>1);
	globalIdx.y -= uint2float((domain_count_grid.y*domain_size_grid.y)>>1);
	
	FLOATd bbox_ll, bbox_lr, bbox_ul, bbox_ur;

	__set_bbbox( &bbox_ll, globalIdx.x-kernel_width/2.0f, globalIdx.y-kernel_width/2.0f, get_last_dim(globalIdx) );
	__set_bbbox( &bbox_lr, globalIdx.x+domain_size_grid.x-1+kernel_width/2.0f, globalIdx.y-kernel_width/2.0f, get_last_dim(globalIdx) );
	__set_bbbox( &bbox_ul, globalIdx.x-kernel_width/2.0f, globalIdx.y+domain_size_grid.y-1+kernel_width/2.0f, get_last_dim(globalIdx) );
	__set_bbbox( &bbox_ur, globalIdx.x+domain_size_grid.x-1+kernel_width/2.0f, globalIdx.y+domain_size_grid.y-1+kernel_width/2.0f,	get_last_dim(globalIdx) );

//	printf("\nbbox_ul: %f %f %f", bbox_ul.x, bbox_ul.y, bbox_ul.z );
//	printf("\nbbox_lr: %f %f %f", bbox_lr.x, bbox_lr.y, bbox_lr.z );

	unsigned char scenario;
	float theta_min, theta_max;

	// Take special care if any of the bbox corners is exactly the origin (atan2f undefined then)
	if( (bbox_ul.x == 0.0f && bbox_ul.y == 0.0f) || (bbox_ur.x == 0.0f && bbox_ur.y == 0.0f) || (bbox_ll.x == 0.0f && bbox_ll.y == 0.0f) || (bbox_lr.x == 0.0f && bbox_lr.y == 0.0f) )
	{
		// Treat as if the origin is fully inside the bbox. All projections are hit.
		scenario = 1; 
	}
	else{

		// Find angles for bounding box corners
		float _theta1 = atan2f(bbox_ul.y,bbox_ul.x);
		float _theta2 = atan2f(bbox_ur.y,bbox_ur.x);
		float _theta3 = atan2f(bbox_ll.y,bbox_ll.x);
		float _theta4 = atan2f(bbox_lr.y,bbox_lr.x);

		//	printf("\ntheta: %f %f %f %f", _theta1, _theta2, _theta3, _theta4 );

		// Find minimum and maximum angles
		float _theta_min1, _theta_min2, _theta_max1, _theta_max2;

		_theta_min1 = min( _theta1, _theta2 );
		_theta_min2 = min( _theta3, _theta4 );
		_theta_max1 = max( _theta1, _theta2 );
		_theta_max2 = max( _theta3, _theta4 );

		theta_min = min( _theta_min1, _theta_min2 );
		theta_max = max( _theta_max1, _theta_max2 );

		if( ((theta_min*theta_max)>=0) || (theta_min>-CUDART_PIO2_F)){

			// Normal case, no discontinuity 
			scenario = 0;
		}
		else if(bbox_ur.x>=0){ 

			// bbox crosses both 'x' and 'y' axis, is thus "around" the origin
			scenario = 1;
		}
		else{

			// bbox crosses 'x' axis with all negative 'x' values. Thus 'theta' has a discontinuity.
			scenario = 2;

			if(_theta1<0)
				_theta1 += (2.0f*CUDART_PI_F);
			if(_theta2<0)
				_theta2 += (2.0f*CUDART_PI_F);
			if(_theta3<0)
				_theta3 += (2.0f*CUDART_PI_F);
			if(_theta4<0)
				_theta4 += (2.0f*CUDART_PI_F);

			_theta_min1 = min( _theta1, _theta2 );
			_theta_min2 = min( _theta3, _theta4 );
			_theta_max1 = max( _theta1, _theta2 );
			_theta_max2 = max( _theta3, _theta4 );

			theta_min = min( _theta_min1, _theta_min2 );
			theta_max = max( _theta_max1, _theta_max2 );
		}
	}

	unsigned int cur_projection = get_last_dim(fixed_dims)*get_last_dim(_globalIdx)*num_projections_per_frame;
	unsigned int num_strips = 0;

	// Test all projections for intersection with bounding box
	while( cur_projection < num_projections_per_frame*(get_last_dim(fixed_dims)*get_last_dim(_globalIdx)+1) ){

		float theta;

		// Calculate angle for current projection [-PI;PI[ (better not accumulate?)

		switch(TYPE){
			case 1:
				{
					const float angle_step = CUDART_PI_F/((sqrtf(5.0f) + 1.0f)/2.0f); // Golden angle step size
					theta = fmodf( (uint2float(cur_projection)+__angular_offset)*angle_step*__gc_factor, 2.0f*CUDART_PI_F );
				}
				break;

			case 2:
				{
					const float cur_proj = uint2float(cur_projection);
					const float frame = floorf((cur_proj+0.5f)*__one_over_num_projections_per_frame);
					const float cur_proj_local = cur_proj - frame*__num_projections_per_frame;
					const float rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
					theta = fmodf( cur_proj_local*CUDART_PI_F*__one_over_num_projections_per_frame+rotation*__interframe_rotation, 2.0f*CUDART_PI_F );
				}
				break;

			default:
				break;
		}

		if( theta>CUDART_PI_F )
			theta -= (2.0f*CUDART_PI_F);

		// Check if theta is inside bounding box (plus mirrored box)
		switch( scenario ){

				case 0:

					// Normal case, no discontinuity 
					if( (theta>theta_min && theta<theta_max) || (theta>(theta_min+CUDART_PI_F) && theta<(theta_max+CUDART_PI_F)) || (theta>(theta_min-CUDART_PI_F) && theta<(theta_max-CUDART_PI_F)) )
						num_strips++;
					break;

				case 1:

					// bbox crosses both 'x' and 'y' axis, is thus "around" the origin
					num_strips++;
					break;

				case 2:

					// bbox crosses 'x' axis with all negative 'x' values. Thus 'theta' has a discontinuity.
					if(theta<0)
						theta += (2.0f*CUDART_PI_F);

					// Recheck if theta is inside bounding box (plus mirrored box). Theta in [0;2*PI[
					if( (theta>theta_min && theta<theta_max) || (theta>(theta_min+CUDART_PI_F) && theta<(theta_max+CUDART_PI_F)) || (theta>(theta_min-CUDART_PI_F) && theta<(theta_max-CUDART_PI_F)) )
						num_strips++;
					break;

				default:
					break;
		}

		cur_projection++;	
	}
	stripCounts[globalThreadId] = num_strips;
	stripCountPairs[globalThreadId] = make_uint2(num_strips,globalThreadId);
}


__global__ void
preprocess_radial_threshold_kernel( unsigned int *stripCounts, unsigned int *threadCounts, unsigned int num_elements )
{
	const unsigned int idx = (blockIdx.x*blockDim.x+threadIdx.x);

	if( idx < num_elements ){
		threadCounts[idx] = (stripCounts[idx]>0);
	}
}

__global__ void
preprocess_radial_compute_domains_map_kernel( unsigned int *threadCounts, unsigned int *threadOffsets, unsigned int *domainsMap, unsigned int num_elements )
{
	const unsigned int idx = (blockIdx.x*blockDim.x+threadIdx.x);

	if( idx < num_elements ){
		if( threadCounts[idx] ){
			domainsMap[threadOffsets[idx]] = idx;
		}
	}
}


template< class UINTd, class FLOATd, char TYPE > void
mr_recon::NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_radial_allocate( unsigned int *num_strips, unsigned int *stripOffsets )
{
	// Use cuDPP to compact the memory required by the NFFT^H

	unsigned int *stripCounts;
	cudaMalloc( (void**) &stripCounts, prod(domain_count_grid)*sizeof(unsigned int) );
	
	if( domainsMapDevPtr_NFFT_H )
		cudaFree(domainsMapDevPtr_NFFT_H);
	cudaMalloc( (void**) &domainsMapDevPtr_NFFT_H, prod(domain_count_grid)*sizeof(uint2) );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_radial_allocate':%s\n", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	// Set kernel constants
	_preproc_NFFT_set_constants( this );

	const unsigned int block_size = 192;
	dim3 blockDim( block_size, 1, 1 );
	dim3 gridDim( (unsigned int) ceil((double)(prod(domain_count_grid)/get_last_dim(domain_count_grid))/(double)block_size), get_last_dim(domain_count_grid) );

    // calculate number of strips and threads needed per domain
	preprocess_radial_allocate_kernel<UINTd,FLOATd,TYPE><<<gridDim, blockDim>>>( matrix_size_os, domain_count_grid, domain_size_grid, fixed_dims, W, stripCounts, domainsMapDevPtr_NFFT_H );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error detected in 'preprocess_radial_allocate_kernel': %s. Quitting.\n", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	int dev;
	cudaGetDevice(&dev);

	static bool _make_plan[] = {true, true, true, true, true, true, true, true};
	static unsigned int _prod_domain_count_grid[] = {prod(domain_count_grid),prod(domain_count_grid),prod(domain_count_grid),prod(domain_count_grid),prod(domain_count_grid),prod(domain_count_grid),prod(domain_count_grid),prod(domain_count_grid)};

	if( _make_plan[dev] || (_prod_domain_count_grid[dev] != prod(domain_count_grid)) ){
		// Scan stripcount array
		CUDPPConfiguration config;
		config.algorithm = CUDPP_SCAN;
		config.datatype = CUDPP_UINT;
		config.op = CUDPP_ADD;
		config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;
		cudppPlan(&scanplan, config, prod(domain_count_grid), 1, 0);

		_prod_domain_count_grid[dev] = prod(domain_count_grid);
		_make_plan[dev] = false;
	}
	cudppScan(scanplan, stripOffsets, stripCounts, prod(domain_count_grid));

	// readback total number of strips
	unsigned int lastElement, lastScanElement;
	cudaMemcpy((void *) &lastElement, 
		(void *) (stripCounts + prod(domain_count_grid)-1), 
		sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy((void *) &lastScanElement, 
		(void *) (stripOffsets + prod(domain_count_grid)-1), 
		sizeof(unsigned int), cudaMemcpyDeviceToHost);
	*num_strips = lastElement + lastScanElement;

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error detected in 'preprocess_radial_allocate': %s. Quitting.\n", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	cudaFree( stripCounts );

	// Radix-sort domains map
	uint2 *__tmp;
	cudaMalloc( (void**) &__tmp, prod(domain_count_grid)*sizeof(uint2) );
	RadixSort((KeyValuePair*)domainsMapDevPtr_NFFT_H, (KeyValuePair*) __tmp, prod(domain_count_grid), 32); 
	cudaFree(__tmp);

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error detected in 'preprocess_radial_allocate': %s. Quitting.\n", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}
}


template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_radial( unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	successfully_preprocessed = false;

	// Some sanity checks

	if( sizeof(UINTd)==sizeof(uint2) && projections_per_frame != num_projections ){
		printf("\nError: 'NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_radial': For 2D reconstructions the projections must reflect a single frame.\n");
		return 0x0;
	}

	if( sizeof(UINTd)<sizeof(uint3) && sum(fixed_dims) ){
		printf("\nERROR: NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_radial : There can be no fixed dimensions in the 2D case.\n");
		return false;				
	}

	if( sizeof(UINTd)==sizeof(uint3) && !get_last_dim(fixed_dims) ){
		printf("\nERROR: NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_radial : For now, the last dimension must be fixed in the 3D case.\n");
		return false;				
	}

	if( sum(fixed_dims)>1 || (sum(fixed_dims)==1 && !get_last_dim(fixed_dims)) ){
		printf("\nERROR: NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_radial : Only one fixed dimension allowed. It must be the last dimension.\n");
		return false;				
	}
/*
	if( prod(domain_count_grid)%preprocess_radial_NFFT_H_block_dim_x ){
		printf("\nERROR: NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_radial : Domain_count.x is not a multiple of the block size.\n");
		return false;
	}
*/
	// Set number of samples
	number_of_samples = num_projections*samples_per_projection;

	// Update plan
	this->samples_per_projection = samples_per_projection;
	this->projections_per_frame = projections_per_frame;
	this->angular_offset = angular_offset;
	this->frames_per_rotation = frames_per_rotation;
	this->interframe_rotation = CUDART_PI_F/(float)(frames_per_rotation*projections_per_frame);
	this->gc_factor = gc_factor;
	this->total_projections_f = (float)num_projections;

	if( domainsMapDevPtr_NFFT_H ){
		cudaFree( domainsMapDevPtr_NFFT_H );
		domainsMapDevPtr_NFFT_H = 0x0;
	}
	if( stripsMapDevPtr_NFFT_H ){
		cudaFree( stripsMapDevPtr_NFFT_H );
		stripsMapDevPtr_NFFT_H = 0x0;
	}
	if( stripsDevPtr_NFFT_H ){
		cudaFree( stripsDevPtr_NFFT_H );
		stripsDevPtr_NFFT_H = 0x0;
	}

	unsigned int *stripOffsets;
	cudaMalloc( (void**) &stripOffsets, prod(domain_count_grid)*sizeof(unsigned int) );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_radial_trajectories_NFFT': %s. Quitting.", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	dim3 blockDim(preprocess_radial_NFFT_H_block_dim_x, 1,1); // Has to be one-dimensional. Arrays in shared memory depend on this!
	dim3 gridDim((unsigned int)ceilf((float)(prod(domain_count_grid)/get_last_dim(domain_count_grid))/(float)blockDim.x), get_last_dim(domain_count_grid) );

	// Compute strip offset map
	unsigned int num_strips/*, num_threads*/;
	preprocess_radial_allocate( &num_strips, stripOffsets );

	// At the center every projection makes a strip. This is the worst case.
	this->number_of_strips_NFFT_H = num_strips;

	cudaMalloc( (void**) &stripsMapDevPtr_NFFT_H, prod(domain_count_grid)*sizeof(uint2) );
	cudaMalloc( (void**) &stripsDevPtr_NFFT_H, num_strips*sizeof(uint2) );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_radial_trajectories_NFFT': %s. Quitting.\n", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	// Set kernel constants
	_preproc_NFFT_set_constants( this );

	// Preprocess
	preprocess_radial_trajectories_NFFT_H_kernel<UINTd,FLOATd,TYPE><<<gridDim, blockDim>>>( matrix_size_os, domain_count_grid, domain_size_grid, fixed_dims, W, stripOffsets, stripsMapDevPtr_NFFT_H, stripsDevPtr_NFFT_H );

	cudaFree( stripOffsets );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'preprocess_radial_trajectories_NFFT_H': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	successfully_preprocessed = true;
		
	return true;
}

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_iteration_plan< UINTd, FLOATd, TYPE >::preprocess_radial( unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	this->successfully_preprocessed = false;

	// Free device memory before re-allocating. 
	// Do this first do reduce the likelyhood of memory shortage...
	
	if( stripsMapDevPtr_NFFT ){
		cudaFree( stripsMapDevPtr_NFFT );
		stripsMapDevPtr_NFFT = 0x0;
	}
	if( stripsDirDevPtr_NFFT ){
		cudaFree( stripsDirDevPtr_NFFT );
		stripsDirDevPtr_NFFT = 0x0;
	}
	if( stripOriginsDevPtr_NFFT ){
		cudaFree( stripOriginsDevPtr_NFFT );
		stripOriginsDevPtr_NFFT = 0x0;
	}
	if( stripLengthsDevPtr_NFFT ){
		cudaFree( stripLengthsDevPtr_NFFT );
		stripLengthsDevPtr_NFFT = 0x0;
	}
	if( domainsMapDevPtr_NFFT_H ){
		cudaFree( domainsMapDevPtr_NFFT_H );
		domainsMapDevPtr_NFFT_H = 0x0;
	}
	if( stripsMapDevPtr_NFFT_H ){
		cudaFree( stripsMapDevPtr_NFFT_H );
		stripsMapDevPtr_NFFT_H = 0x0;
	}
	if( stripsDevPtr_NFFT_H ){
		cudaFree( stripsDevPtr_NFFT_H );
		stripsDevPtr_NFFT_H = 0x0;
	}

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'NFFT_iteration_plan::preprocess_radial': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	NFFT_plan< UINTd, FLOATd, TYPE > *pre_NFFT = new NFFT_plan< UINTd, FLOATd, TYPE >( matrix_size, matrix_size_os, fixed_dims, domain_size_samples, domain_size_coils, W );
	NFFT_H_plan< UINTd, FLOATd, TYPE > *pre_NFFT_H = new NFFT_H_plan< UINTd, FLOATd, TYPE >( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_coils, W );

	pre_NFFT->preprocess_radial( num_projections, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );
	pre_NFFT_H->preprocess_radial( num_projections, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );

	if( !pre_NFFT->successfully_preprocessed || !pre_NFFT_H->successfully_preprocessed ){
		NFFT_cleanup(&pre_NFFT);
		NFFT_cleanup(&pre_NFFT_H);
		return false;
	}

	assert( pre_NFFT->samples_per_projection == pre_NFFT_H->samples_per_projection );
	assert( pre_NFFT->projections_per_frame == pre_NFFT_H->projections_per_frame );
	assert( pre_NFFT->angular_offset == pre_NFFT_H->angular_offset );
	assert( pre_NFFT->frames_per_rotation == pre_NFFT_H->frames_per_rotation );
	assert( pre_NFFT->interframe_rotation == pre_NFFT_H->interframe_rotation );
	assert( pre_NFFT->gc_factor == pre_NFFT_H->gc_factor );
	assert( pre_NFFT->total_projections_f == pre_NFFT_H->total_projections_f );

	this->samples_per_projection = pre_NFFT->samples_per_projection;
	this->projections_per_frame = pre_NFFT->projections_per_frame;
	this->angular_offset = pre_NFFT->angular_offset;
	this->frames_per_rotation = pre_NFFT->frames_per_rotation;
	this->interframe_rotation = pre_NFFT->interframe_rotation;
	this->gc_factor = pre_NFFT->gc_factor;
	this->total_projections_f = pre_NFFT->total_projections_f;

	assert( pre_NFFT->number_of_samples == pre_NFFT_H->number_of_samples );
	this->number_of_samples = pre_NFFT->number_of_samples;

	this->domain_count_samples = pre_NFFT->domain_count_samples;

	this->max_strips_per_domain_NFFT = pre_NFFT->max_strips_per_domain_NFFT;

	this->number_of_strips_NFFT = pre_NFFT->number_of_strips_NFFT;
//	this->number_of_threads_NFFT_H = pre_NFFT_H->number_of_threads_NFFT_H;
	this->number_of_strips_NFFT_H = pre_NFFT_H->number_of_strips_NFFT_H;

	this->stripsMapDevPtr_NFFT = pre_NFFT->stripsMapDevPtr_NFFT;
	this->stripsDirDevPtr_NFFT = pre_NFFT->stripsDirDevPtr_NFFT;
	this->stripOriginsDevPtr_NFFT = pre_NFFT->stripOriginsDevPtr_NFFT;
	this->stripLengthsDevPtr_NFFT = pre_NFFT->stripLengthsDevPtr_NFFT;

	this->domainsMapDevPtr_NFFT_H = pre_NFFT_H->domainsMapDevPtr_NFFT_H;
	this->stripsMapDevPtr_NFFT_H = pre_NFFT_H->stripsMapDevPtr_NFFT_H;
	this->stripsDevPtr_NFFT_H = pre_NFFT_H->stripsDevPtr_NFFT_H;

	delete pre_NFFT;
	delete pre_NFFT_H;

	this->successfully_preprocessed = true;

	return true;
}


// Find the (eight) neighbors to a given radial sample index


template <char TYPE> __inline__ __device__ float2
compute_radial_neighbors( unsigned int sampleIdx, int index_on_projection, float alpha, unsigned int bias, float2 *p1, float2 *p2, float2 *p3, float2 *p4, float2 *p5, float2 *p6, float2 *p7, float2 *p8  )
{
	// This code assumes that all eight neighbors exist. 
	// I.e. don't invoke on an index that is the first or last sample of a projection.

	float cos_angle, sin_angle;

	const unsigned int num_samples_per_projection = float2uint(__num_samples_per_projection);
	const unsigned int num_projections_per_frame = float2uint(__num_projections_per_frame);

	// The projections indices we will investigate
	const unsigned int current_projection = sampleIdx / num_samples_per_projection;
	const unsigned int current_frame = current_projection / num_projections_per_frame;
	
	// The sample positions (scales) can be either of the _local_ indices 'index_on_projection' or 'samples_per_projection'-'index_on_projection'
	// Beware of "skewness" around the origin, i.e. +1 sample on one side
	const float ctr_scale		= alpha*(index_on_projection*__one_over_radial_oversampling-bias);
	const float ctr_scale_inv	= alpha*((__num_samples_per_projection-index_on_projection)*__one_over_radial_oversampling-bias);
	const float prev_scale		= alpha*((index_on_projection-1)*__one_over_radial_oversampling-bias);
	const float prev_scale_inv	= alpha*((__num_samples_per_projection-(index_on_projection-1))*__one_over_radial_oversampling-bias);
	const float next_scale		= alpha*((index_on_projection+1)*__one_over_radial_oversampling-bias);
	const float next_scale_inv	= alpha*((__num_samples_per_projection-(index_on_projection+1))*__one_over_radial_oversampling-bias);

	// Unit circle position for current projection
	switch(TYPE){
	
		case 1:
			{
				const float angle_step = CUDART_PI_F/((sqrtf(5.0f) + 1.0f)/2.0f); // Golden ratio
				__sincosf( (uint2float(current_projection)+__angular_offset)*angle_step*__gc_factor, &sin_angle, &cos_angle );
			}
			break;

		case 2:
			{
				const float cur_proj = uint2float(current_projection);
				const float frame = floorf((cur_proj+0.5f)*__one_over_num_projections_per_frame);
				const float cur_proj_local = cur_proj - frame*__num_projections_per_frame;
				const float rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
				__sincosf( cur_proj_local*CUDART_PI_F*__one_over_num_projections_per_frame+rotation*__interframe_rotation, &sin_angle, &cos_angle );
			}
			break;

		default: 
			break;
	}

	// Find the normal to the current projection direction
	float2 normal = make_float2( -sin_angle, cos_angle );

	// The position of the sampleIdx itself
	const float2 sample_pos = make_float2( ctr_scale*cos_angle, ctr_scale*sin_angle );

	// The positions of the previous and next sample
	*p1 = make_float2( prev_scale*cos_angle, prev_scale*sin_angle );
	*p2 = make_float2( next_scale*cos_angle, next_scale*sin_angle );

	// Initialize remaining points;
	*p3 = *p4 = *p5 = *p6 = *p7 = *p8 = make_float2( 999999999.9f, 999999999.9f ); // Far away...

	// Run through all projections to find the closests neighbors
	
	for( unsigned int i=current_frame*num_projections_per_frame; i<(current_frame+1)*__num_projections_per_frame; i++ ){

		if( i == current_projection )
			continue;

		// Unit circle position projection 'i'
		switch(TYPE){

		case 1:
			{
				const float angle_step = CUDART_PI_F/((sqrtf(5.0f) + 1.0f)/2.0f); // Golden ratio
				__sincosf( (uint2float(i)+__angular_offset)*angle_step*__gc_factor, &sin_angle, &cos_angle );
			}
			break;

		case 2:
			{
				const float cur_proj = uint2float(i);
				const float frame = floorf((cur_proj+0.5f)*__one_over_num_projections_per_frame);
				const float cur_proj_local = cur_proj - frame*__num_projections_per_frame;
				const float rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
				__sincosf( cur_proj_local*CUDART_PI_F*__one_over_num_projections_per_frame+rotation*__interframe_rotation, &sin_angle, &cos_angle );
			}
			break;

		default: 
			break;
		}

		// Determine sample positions on projection
		float2 prev_pos_1 = make_float2( prev_scale*cos_angle,		prev_scale*sin_angle );
		float2 prev_pos_2 = make_float2( prev_scale_inv*cos_angle,	prev_scale_inv*sin_angle );
		float2 ctr_pos_1  = make_float2( ctr_scale*cos_angle,		ctr_scale*sin_angle );
		float2 ctr_pos_2  = make_float2( ctr_scale_inv*cos_angle,	ctr_scale_inv*sin_angle );
		float2 next_pos_1 = make_float2( next_scale*cos_angle,		next_scale*sin_angle );
		float2 next_pos_2 = make_float2( next_scale_inv*cos_angle,	next_scale_inv*sin_angle );

		// The dot product is used to ensure we find a neighbor on each side
		if( dot(ctr_pos_1-sample_pos, normal) > 0 ){

			if( squared_length(ctr_pos_1-sample_pos) < squared_length(*p4-sample_pos) ){
				*p3 = prev_pos_1;
				*p4 = ctr_pos_1;
				*p5 = next_pos_1;
			}
		}
		else{

			if( squared_length(ctr_pos_1-sample_pos) < squared_length(*p7-sample_pos) ){
				*p6 = prev_pos_1;
				*p7 = ctr_pos_1;
				*p8 = next_pos_1;
			}
		}
	
		// The dot product is used to ensure we find a neighbor on each side
		if( dot(ctr_pos_2-sample_pos, normal) > 0 ){

			if( squared_length(ctr_pos_2-sample_pos) < squared_length(*p4-sample_pos) ){
				*p3 = prev_pos_2;
				*p4 = ctr_pos_2;
				*p5 = next_pos_2;
			}
		}
		else{

			if( squared_length(ctr_pos_2-sample_pos) < squared_length(*p7-sample_pos) ){
				*p6 = prev_pos_2;
				*p7 = ctr_pos_2;
				*p8 = next_pos_2;
			}
		}
	}
	
	return sample_pos;
}


// Compute density compensation weights, radial trajectory specialization

template <char TYPE> __global__ void 
compute_dcw_radial_kernel( float alpha, unsigned int number_of_samples, unsigned int bias, float *weights_DevPtr )
{
	// Computes in the special (easy) case of radial trajectories

	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_samples_per_projection = float2uint(__num_samples_per_projection);
	const unsigned int num_projections_per_frame = float2uint(__num_projections_per_frame);

	if( index < number_of_samples ){

		float weight;
		const int index_on_projection = index%num_samples_per_projection;

		if( index_on_projection == (num_samples_per_projection>>1) ){
			// Special case - center of k-space
			const float radius = (alpha*__one_over_radial_oversampling)/2.0f;
			const float area = radius*radius*CUDART_PI_F;
			weight = area/__num_projections_per_frame;
		}
//		else if( (index_on_projection == 0) || (index_on_projection == (samples_per_projection-1)) ){
//			// Special case - outermost sample in projection/k-space
//			weight = 256;
//		}
		else{

			// General case - all neighbors exist
			
			// Compute sample positions for the current sample and all neighbors
			// The ordering of p1..p8 in the call below follows the edge of the "Voronoi polygon"

			float2 sample_pos;
			float2 p1, p2, p3, p4, p5, p6, p7, p8;

			sample_pos = compute_radial_neighbors<TYPE>( index, index_on_projection, alpha, bias, &p1, &p5, &p2, &p3, &p4, &p8, &p7, &p6 );

			// Find midpoints of lines from sample_pos to all other points.
			p1 = 0.5f*(sample_pos+p1);
			p2 = 0.5f*(sample_pos+p2);
			p3 = 0.5f*(sample_pos+p3);
			p4 = 0.5f*(sample_pos+p4);
			p5 = 0.5f*(sample_pos+p5);
			p6 = 0.5f*(sample_pos+p6);
			p7 = 0.5f*(sample_pos+p7);
			p8 = 0.5f*(sample_pos+p8);

			// The weight is determined by the area of the polygon (http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/)
			weight = 0.5f*((p1.x*p2.y-p2.x*p1.y)+(p2.x*p3.y-p3.x*p2.y)+(p3.x*p4.y-p4.x*p3.y)+(p4.x*p5.y-p5.x*p4.y)+(p5.x*p6.y-p6.x*p5.y)+(p6.x*p7.y-p7.x*p6.y)+(p7.x*p8.y-p8.x*p7.y)+(p8.x*p1.y-p1.x*p8.y));			
			if( weight<0 ) weight *= -1.0f;
		}

		weights_DevPtr[index] = weight;
	}
}

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN  > bool
preprocess_radial_NFFT( PLAN<UINTd, FLOATd, TYPE> *plan, unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	return plan->preprocess_radial( num_projections, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );
}

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > bool
compute_dcw_radial( PLAN<UINTd, FLOATd, TYPE> *plan )
{
	if( plan->projections_per_frame < 4 ){
		printf("\n'Warning : compute_dcw_radial' : Please use at least four projections per frame for proper computation of weights. Returning.");
		return false;
	}

	// Allocate device memory
	if( plan->weights_DevPtr )
		cudaFree( plan->weights_DevPtr );
	cudaMalloc( (void**) &plan->weights_DevPtr, plan->number_of_samples*sizeof(float) );

	// Find dimensions of grid/blocks.
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );
	
//	const unsigned int block_size = deviceProp.maxThreadsPerBlock>>2;
	const unsigned int block_size = 192;

	dim3 blockDim( block_size, 1, 1 );
	dim3 gridDim( (unsigned int) ceil((double)plan->number_of_samples/block_size), 1, 1 );

	assert( plan->number_of_samples >= block_size );

	unsigned int bias = plan->matrix_size.x>>1;

	// Invoke kernel
	_preproc_NFFT_set_constants( plan );
	compute_dcw_radial_kernel<TYPE><<<gridDim, blockDim>>>( plan->alpha, plan->number_of_samples, bias, plan->weights_DevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'compute_dcw_radial_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	/*
	//DEBUG
	printf("\nWriting dcw_GPU.raw with %d elements.\n", plan->number_of_samples );
	float* tmp = (float*) calloc( 1, plan->number_of_samples*sizeof(float) );
	cudaMemcpy( tmp, plan->weights_DevPtr, plan->number_of_samples*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("dcw_GPU.raw", "wb");
	fwrite( tmp, plan->number_of_samples, sizeof(float), fout );
	fclose(fout);
	free(tmp);
	*/

	return true;
}

template< char TYPE > float* 
compute_dcw_radial_2d( unsigned int matrix_size, unsigned int matrix_size_os, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	if( projections_per_frame < 4 ){
		printf("\n'Warning : compute_dcw_radial' : Use at least four projections per frame for proper computation of weights. Returning 0x0.");
		return 0x0;
	}

	const unsigned int number_of_samples = samples_per_projection*projections_per_frame;
	const float alpha = (float)matrix_size_os/(float)matrix_size;
	const unsigned int bias = matrix_size>>1;

	// Allocate device memory for dcw
	float *weights_DevPtr;
	cudaMalloc( (void**) &weights_DevPtr, number_of_samples*sizeof(float) );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'compute_dcw_radial': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Find dimensions of grid/blocks.
	const unsigned int block_size = min(192,number_of_samples);
	dim3 blockDim( block_size, 1, 1 );
	dim3 gridDim( (unsigned int) ceil((double)number_of_samples/blockDim.x), 1, 1 );

	// Invoke kernel	
	_preproc_NFFT_set_constants<TYPE>( matrix_size, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );
	compute_dcw_radial_kernel<TYPE><<<gridDim, blockDim>>>( alpha, number_of_samples, bias, weights_DevPtr );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'compute_dcw_radial_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	return weights_DevPtr;
}

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
float* get_dcw( PLAN<UINTd, FLOATd, TYPE> *plan )
{
	return plan->weights_DevPtr;
}

template< char TYPE > float2* 
compute_trajectory_radial_2d( unsigned int matrix_size, unsigned int matrix_size_os, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	const unsigned int number_of_samples = samples_per_projection*projections_per_frame;
	const float alpha = (float)matrix_size_os/(float)matrix_size;
	const unsigned int bias = matrix_size>>1;

	// Find dimensions of grid/blocks.
	const unsigned int block_size = min(256,number_of_samples);
	dim3 dimBlock( block_size, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)number_of_samples/dimBlock.x), 1, 1 );

	// Invoke kernel
	float2 *co = compute_radial_sample_positions<uint2,float2,TYPE>( number_of_samples, make_uint2(matrix_size,matrix_size), make_uint2(matrix_size_os,matrix_size_os), alpha, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor, true, true );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'compute_trajectory_radial_2d': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	return co;
}


/*
	Noise scaling (decorrelation).
	Output a dataset matching the data samples size of noise magnitudes squared outside shutter and 0 inside.
*/

// NOTICE: There is a separat version of this code in preprocess_kt.cu
//		   Any changes made here is likely necessary in the general case also.

template <char TYPE > __global__ void
noise_decorrelate_radial_kernel( unsigned int num_points, unsigned int num_coils, float shutter_radius_squared, unsigned int bias, cuFloatComplex *in_samples_devPtr, float *out_magnitudes_squared, float *out_count  )
{
	const unsigned int sample_no = (blockIdx.x*blockDim.x+threadIdx.x);
	const unsigned int num_samples_per_projection = float2uint(__num_samples_per_projection);

	if( sample_no >= num_points )
		return;

	/*
		Compute sample position. 
	*/

	float cos_angle, sin_angle;

	// The projection indices we will investigate
	const unsigned int index_on_projection = sample_no%num_samples_per_projection;
	const unsigned int current_projection = sample_no/num_samples_per_projection;
	const float scale = index_on_projection*__one_over_radial_oversampling-bias;

	switch(TYPE){

	case 1:
		{
			// Unit circle position for current projection
			const float angle_step = CUDART_PI_F/((sqrtf(5.0f) + 1.0f)/2.0f); // Golden ratio
			__sincosf( (uint2float(current_projection)+__angular_offset)*angle_step*__gc_factor, &sin_angle, &cos_angle );
		}
		break;
	
	case 2:
		{
			const float cur_proj = uint2float(current_projection);
			const float frame = floorf((cur_proj+0.5f)*__one_over_num_projections_per_frame);
			const float cur_proj_local = cur_proj - frame*__num_projections_per_frame;
			const float rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
			__sincosf( cur_proj_local*CUDART_PI_F*__one_over_num_projections_per_frame+rotation*__interframe_rotation, &sin_angle, &cos_angle );
		}
		break;

	default:
		break;
	}

	// The position of the sampleIdx itself
	const float2 sample_pos = make_float2( scale*cos_angle, scale*sin_angle );

	float dist_sq = dot(sample_pos, sample_pos);

	// Write to global memory

	if( dist_sq >= shutter_radius_squared ){
		for( unsigned int c=0; c<num_coils; c++ ){
			cuFloatComplex _data = in_samples_devPtr[c*num_points+sample_no];
			out_magnitudes_squared[c*num_points+sample_no] = cuCrealf(cuCmulf(_data, cuConjf(_data)));
		}
		out_count[sample_no] = 1.0f;
	}
	else{
		for( unsigned int c=0; c<num_coils; c++ )
			out_magnitudes_squared[c*num_points+sample_no] = 0.0f;
		out_count[sample_no] = 0.0f;
	}
}

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > bool
noise_decorrelate_radial( PLAN<UINTd, FLOATd, TYPE> *plan, float shutter_radius, cuFloatComplex *samples_DevPtr )
{

	// First find noise magnitudes
	
	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)plan->number_of_samples/blockDim.x), 1, 1 );

	assert(plan->number_of_samples>=blockDim.x);

	float *noise_modulus_sq, *point_count;

	cudaMalloc( (void**) &noise_modulus_sq, plan->number_of_samples*plan->number_of_coils*sizeof(float) );
	cudaMalloc( (void**) &point_count, plan->number_of_samples*sizeof(float) );

	_preproc_NFFT_set_constants( plan );
	noise_decorrelate_radial_kernel<TYPE><<< gridDim, blockDim >>>( plan->number_of_samples, plan->number_of_coils, shutter_radius*shutter_radius, plan->matrix_size.x>>1, samples_DevPtr, noise_modulus_sq, point_count );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'noise_decorrelate_radial_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Then sum for each coil and normalize

	float *noise_variances = (float*) malloc(plan->number_of_coils*sizeof(float));
	float used_points = cublasSasum( plan->number_of_samples, point_count, 1 );

	if( used_points > 0.1f ){
		for( unsigned int c=0; c<plan->number_of_coils; c++ ){
			noise_variances[c] = cublasSasum( plan->number_of_samples, &noise_modulus_sq[c*plan->number_of_samples], 1 );	
			noise_variances[c] /= used_points;
		}
	}
	else{
		printf("\nnoise_scaling: No points found within k-space shutter!\n");
		return false;
	}

	// And finally scale the samples values accordingly

	for( unsigned int c=0; c<plan->number_of_coils; c++){
//		printf("\nNoise scaling with a factor of %f for coil %d", 1.0f/sqrtf(noise_variances[c]), c );
		cublasSscal( 2*plan->number_of_samples, 1.0f/sqrtf(noise_variances[c]), (float*)(&samples_DevPtr[c*plan->number_of_samples]), 1 );
	}

	free(noise_variances);
	cudaFree( noise_modulus_sq );
	cudaFree( point_count );
	
	return true;
}



/*
	Template instantiation
*/

// Only initialize for TYPE 1 and 2 (radial)

template mr_recon::NFFT_plan< uint2, float2, 1 >* preprocess_radial_NFFT( uint2, uint2, uint2, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_H_plan< uint2, float2, 1 >* preprocess_radial_NFFT( uint2, uint2, uint2, uint2, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_iteration_plan< uint2, float2, 1 >* preprocess_radial_NFFT( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template mr_recon::NFFT_plan< uint2, float2, 2 >* preprocess_radial_NFFT( uint2, uint2, uint2, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_H_plan< uint2, float2, 2 >* preprocess_radial_NFFT( uint2, uint2, uint2, uint2, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_iteration_plan< uint2, float2, 2 >* preprocess_radial_NFFT( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template mr_recon::NFFT_plan< uint3, float3, 1 >* preprocess_radial_NFFT( uint3, uint3, uint3, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_H_plan< uint3, float3, 1 >* preprocess_radial_NFFT( uint3, uint3, uint3, uint3, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_iteration_plan< uint3, float3, 1 >* preprocess_radial_NFFT( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template mr_recon::NFFT_plan< uint3, float3, 2 >* preprocess_radial_NFFT( uint3, uint3, uint3, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_H_plan< uint3, float3, 2 >* preprocess_radial_NFFT( uint3, uint3, uint3, uint3, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_iteration_plan< uint3, float3, 2 >* preprocess_radial_NFFT( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template bool preprocess_radial_NFFT( mr_recon::NFFT_plan< uint2, float2, 1>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_H_plan< uint2, float2, 1>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template bool preprocess_radial_NFFT( mr_recon::NFFT_plan< uint2, float2, 2>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_H_plan< uint2, float2, 2>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template bool preprocess_radial_NFFT( mr_recon::NFFT_plan< uint3, float3, 1>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_H_plan< uint3, float3, 1>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template bool preprocess_radial_NFFT( mr_recon::NFFT_plan< uint3, float3, 2>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_H_plan< uint3, float3, 2>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template bool preprocess_radial_NFFT( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template bool compute_dcw_radial( mr_recon::NFFT_H_plan< uint2, float2, 1>* );
template bool compute_dcw_radial( mr_recon::NFFT_iteration_plan< uint2, float2, 1>* );

template bool compute_dcw_radial( mr_recon::NFFT_H_plan< uint2, float2, 2>* );
template bool compute_dcw_radial( mr_recon::NFFT_iteration_plan< uint2, float2, 2>* );

template bool compute_dcw_radial( mr_recon::NFFT_H_plan< uint3, float3, 1>* );
template bool compute_dcw_radial( mr_recon::NFFT_iteration_plan< uint3, float3, 1>* );

template bool compute_dcw_radial( mr_recon::NFFT_H_plan< uint3, float3, 2>* );
template bool compute_dcw_radial( mr_recon::NFFT_iteration_plan< uint3, float3, 2>* );

template float* compute_dcw_radial_2d<1>( unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template float* compute_dcw_radial_2d<2>( unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template float* get_dcw( mr_recon::NFFT_H_plan< uint2, float2, 1>* );
template float* get_dcw( mr_recon::NFFT_iteration_plan< uint2, float2, 1>* );

template float* get_dcw( mr_recon::NFFT_H_plan< uint2, float2, 2>* );
template float* get_dcw( mr_recon::NFFT_iteration_plan< uint2, float2, 2>* );

template float* get_dcw( mr_recon::NFFT_H_plan< uint3, float3, 1>* );
template float* get_dcw( mr_recon::NFFT_iteration_plan< uint3, float3, 1>* );

template float* get_dcw( mr_recon::NFFT_H_plan< uint3, float3, 2>* );
template float* get_dcw( mr_recon::NFFT_iteration_plan< uint3, float3, 2>* );

template float2* compute_trajectory_radial_2d<1>( unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template float2* compute_trajectory_radial_2d<2>( unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template bool noise_decorrelate_radial( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, float, cuFloatComplex* );
template bool noise_decorrelate_radial( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, float, cuFloatComplex* );

template bool noise_decorrelate_radial( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, float, cuFloatComplex* );
template bool noise_decorrelate_radial( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, float, cuFloatComplex* );
