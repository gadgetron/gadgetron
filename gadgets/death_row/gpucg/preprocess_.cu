#include "preprocess.hcu"
#include "preprocess_private.hcu"
#include "uint_util.hcu"
#include "uint_util_device.hcu"
#include "float_util.hcu"
#include "NFFT.hcu"

#include <stdio.h>
#include <iostream>

#include <cublas.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/generate.h>
#include <thrust/pair.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>


#include <assert.h>

using namespace thrust;

template <class FLOATd>
struct compute_num_cells_per_sample
{
	__host__ __device__
		compute_num_cells_per_sample(unsigned int _d, float _half_W) : d(_d), half_W(_half_W) {}

	__host__ __device__
		unsigned int operator()(FLOATd p) const
	{
		unsigned int num_cells = 1;
		for( unsigned int dim=0; dim<d; dim++ ){
				unsigned int upper_limit = (unsigned int)floor((((float*)&p)[dim])+half_W);
			unsigned int lower_limit = (unsigned int)ceil((((float*)&p)[dim])-half_W);
			num_cells *= (upper_limit-lower_limit+1);
		}
		return num_cells;
	}

	unsigned int d;
	float half_W;
};

// convert input trajectory in [-1/2;1/2] to [0;matrix_size_os+matrix_size_wrap]
template<class FLOATd> struct traj_scale
{
    const FLOATd matrix, bias;

    traj_scale(FLOATd m, FLOATd b) : matrix(m), bias(b) {}

    __host__ __device__
        FLOATd operator()(const FLOATd &in) const { 
            return (in*matrix)+bias;
        }
};

__device__ void
output_pairs( unsigned int sample_idx, float2 p, uint2 matrix_size_os, uint2 matrix_size_wrap, float half_W, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{
	unsigned int lower_limit_x = (unsigned int)ceil(p.x-half_W);
	unsigned int lower_limit_y = (unsigned int)ceil(p.y-half_W);
	unsigned int upper_limit_x = (unsigned int)floor(p.x+half_W);
	unsigned int upper_limit_y = (unsigned int)floor(p.y+half_W);

	unsigned int pair_idx = 0;
	unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
	for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
		for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
			tuples_first[write_offset+pair_idx] = co_to_idx(make_uint2(x,y), matrix_size_os+matrix_size_wrap);
			tuples_last[write_offset+pair_idx] = sample_idx;
			pair_idx++;
		}
	}
}

__device__ void
output_pairs( unsigned int sample_idx, float3 p, uint3 matrix_size_os, uint3 matrix_size_wrap, float half_W, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{
	unsigned int lower_limit_x = (unsigned int)ceil(p.x-half_W);
	unsigned int lower_limit_y = (unsigned int)ceil(p.y-half_W);
	unsigned int lower_limit_z = (unsigned int)ceil(p.z-half_W);
	unsigned int upper_limit_x = (unsigned int)floor(p.x+half_W);
	unsigned int upper_limit_y = (unsigned int)floor(p.y+half_W);
	unsigned int upper_limit_z = (unsigned int)floor(p.z+half_W);

	unsigned int pair_idx = 0;
	unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
	for( unsigned int z=lower_limit_z; z<=upper_limit_z; z++ ){
		for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
			for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
				tuples_first[write_offset+pair_idx] = co_to_idx(make_uint3(x,y,z), matrix_size_os+matrix_size_wrap);
				tuples_last[write_offset+pair_idx] = sample_idx;
				pair_idx++;
			}
		}
	}
}

__device__ void
output_pairs( unsigned int sample_idx, float4 p, uint4 matrix_size_os, uint4 matrix_size_wrap, float half_W, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{
	unsigned int lower_limit_x = (unsigned int)ceil(p.x-half_W);
	unsigned int lower_limit_y = (unsigned int)ceil(p.y-half_W);
	unsigned int lower_limit_z = (unsigned int)ceil(p.z-half_W);
	unsigned int lower_limit_w = (unsigned int)ceil(p.w-half_W);
	unsigned int upper_limit_x = (unsigned int)floor(p.x+half_W);
	unsigned int upper_limit_y = (unsigned int)floor(p.y+half_W);
	unsigned int upper_limit_z = (unsigned int)floor(p.z+half_W);
	unsigned int upper_limit_w = (unsigned int)floor(p.w+half_W);

	unsigned int pair_idx = 0;
	unsigned int write_offset = (sample_idx==0) ? 0 : write_offsets[sample_idx-1];
	for( unsigned int w=lower_limit_w; w<=upper_limit_w; w++ ){
		for( unsigned int z=lower_limit_z; z<=upper_limit_z; z++ ){
			for( unsigned int y=lower_limit_y; y<=upper_limit_y; y++ ){
				for( unsigned int x=lower_limit_x; x<=upper_limit_x; x++ ){
					tuples_first[write_offset+pair_idx] = co_to_idx(make_uint4(x,y,z,w), matrix_size_os+matrix_size_wrap);
					tuples_last[write_offset+pair_idx] = sample_idx;
					pair_idx++;
				}
			}
		}
	}
}

template <class UINTd, class FLOATd> __global__ void
write_pairs_kernel( UINTd matrix_size_os, UINTd matrix_size_wrap, unsigned int num_samples, float half_W, FLOATd *traj_positions, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{
	// Get sample idx
	unsigned int sample_idx = blockIdx.x*blockDim.x + threadIdx.x;

	if( sample_idx<num_samples ){

		FLOATd p = traj_positions[sample_idx];
		output_pairs( sample_idx, p, matrix_size_os, matrix_size_wrap, half_W, write_offsets, tuples_first, tuples_last );
	}
};

template <class UINTd, class FLOATd> void 
write_pairs( UINTd matrix_size_os, UINTd matrix_size_wrap, unsigned int num_samples, float W, FLOATd *traj_positions, unsigned int *write_offsets, unsigned int *tuples_first, unsigned int *tuples_last )
{  
	dim3 blockDim(512,1,1);
	dim3 gridDim((int)ceil((double)num_samples/(double)blockDim.x),1,1);

	write_pairs_kernel<UINTd, FLOATd><<< gridDim, blockDim >>>( matrix_size_os, matrix_size_wrap, num_samples, W/2.0f, traj_positions, write_offsets, tuples_first, tuples_last );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'write_pairs': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Generic preprocessing
template< class UINTd, class FLOATd > mr_recon::NFFT_H_plan<UINTd, FLOATd, 0 >*
preprocess_generic_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_coils, float W, unsigned int num_samples, FLOATd *traj )
{
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	mr_recon::NFFT_H_plan< UINTd, FLOATd, 0 > *tmp = new mr_recon::NFFT_H_plan< UINTd, FLOATd, 0 >( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_coils, W );

	tmp->preprocess_generic( num_samples, traj );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}

template< class UINTd, class FLOATd, template< class, class, char > class PLAN  > bool
preprocess_generic_NFFT( PLAN<UINTd, FLOATd, 0> *plan, unsigned int num_samples, FLOATd *traj )
{
	return plan->preprocess_generic( num_samples, traj );
}

template< class UINTd, class FLOATd > mr_recon::NFFT_iteration_plan<UINTd, FLOATd, 0 >*
preprocess_generic_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, unsigned int num_samples, FLOATd *traj )
{
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	mr_recon::NFFT_iteration_plan<UINTd, FLOATd, 0 > *tmp = new mr_recon::NFFT_iteration_plan<UINTd, FLOATd, 0 >( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_samples, domain_size_coils, W );

	tmp->preprocess_generic( num_samples, traj );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}

void make_bias_vec(float bias, float2 *ret)
{
	ret->x = ret->y = bias;
}
void make_bias_vec(float bias, float3 *ret)
{
	ret->x = ret->y = ret->z = bias;
}
void make_bias_vec(float bias, float4 *ret)
{
	ret->x = ret->y = ret->z = ret->w = bias;
}

inline __host__ float2 uintd_to_floatd_h( uint2 a )
{
	return make_float2( a.x, a.y );
}

inline __host__ float3 uintd_to_floatd_h( uint3 a )
{
	return make_float3( a.x, a.y, a.z );
}

inline __host__ float4 uintd_to_floatd_h( uint4 a )
{
	return make_float4( a.x, a.y, a.z, a.w );
}

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_plan< UINTd, FLOATd, TYPE >::preprocess_generic( unsigned int num_samples, FLOATd *traj )
{
	// TODO...
	return false;
}

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess_generic( unsigned int num_samples, FLOATd *traj )
{
	successfully_preprocessed = false;

	number_of_samples = num_samples;

	// Some sanity checks

	if( sizeof(UINTd)<sizeof(uint3) && sum(fixed_dims) ){
		printf("\nERROR: NFFT_H_plan< UINTd, FLOATd, 0 >::preprocess_generic : There can be no fixed dimensions in the 2D case.\n");
		return false;				
	}

	if( sum(fixed_dims)>1 || (sum(fixed_dims)==1 && !get_last_dim(fixed_dims)) ){
		printf("\nERROR: NFFT_H_plan< UINTd, FLOATd, 0 >::preprocess_generic : Only one fixed dimension allowed. It must be the last dimension.\n");
		return false;				
	}

	// Nullify radial only entries
	this->samples_per_projection = 0;
	this->projections_per_frame = 0;
	this->angular_offset = 0;
	this->frames_per_rotation = 0;
	this->interframe_rotation = 0;
	this->gc_factor = 0;
	this->total_projections_f = 0;
	this->number_of_strips_NFFT_H = 0;

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

	if(	traj_positions ){
		delete traj_positions;
		traj_positions = 0x0;
	}

	if(	tuples_last ){
		delete tuples_last;
		tuples_last = 0x0;
	}

	if(	bucket_begin ){
		delete bucket_begin;
		bucket_begin = 0x0;
	}

	if(	bucket_end ){
		delete bucket_end;
		bucket_end = 0x0;
	}

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_generic_NFFT': %s. Quitting.", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	// Preprocess
	// make Thrust device vector of trajectory and samples
	device_vector<FLOATd> traj_positions_in( device_pointer_cast<FLOATd>(traj), device_pointer_cast<FLOATd>(traj+num_samples) );
	traj_positions = new device_vector<FLOATd>( num_samples );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_generic_NFFT': %s. Quitting.", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	// convert input trajectory in [-1/2;1/2] to [0;matrix_size_os]
	transform( traj_positions_in.begin(), traj_positions_in.end(), traj_positions->begin(), traj_scale<FLOATd>(uintd_to_floatd_h(matrix_size_os), uintd_to_floatd_h((matrix_size_os+matrix_size_wrap)>>1)) );
	
	// allocate storage for and compute temporary prefix-sum variable (#cells influenced per sample)
	device_vector<unsigned int> c_p_s(num_samples);
	device_vector<unsigned int> c_p_s_ps(num_samples);

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_generic_NFFT': %s. Quitting.", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	transform(traj_positions->begin(), traj_positions->end(), c_p_s.begin(), compute_num_cells_per_sample<FLOATd>(d,W/2.0f));
	thrust::plus<unsigned int> binary_op;
	inclusive_scan( c_p_s.begin(), c_p_s.end(), c_p_s_ps.begin(), binary_op ); // prefix sum

	// Build the vector of (grid_idx, sample_idx) tuples. Actually kept in two seperate vectors.
	unsigned int num_pairs = c_p_s_ps.back();
	thrust::device_vector<unsigned int> tuples_first = device_vector<unsigned int>(num_pairs);
	tuples_last = new device_vector<unsigned int>(num_pairs);

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_generic_NFFT': %s. Quitting.", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	// Fill tuple vector
	write_pairs<UINTd, FLOATd>( matrix_size_os, matrix_size_wrap, number_of_samples, W, raw_pointer_cast(&(*traj_positions)[0]), raw_pointer_cast(&c_p_s_ps[0]), raw_pointer_cast(&tuples_first[0]), raw_pointer_cast(&(*tuples_last)[0]) );

	// Sort by grid indices
	sort_by_key(tuples_first.begin(), tuples_first.end(), tuples_last->begin() );

	// each bucket_begin[i] indexes the first element of bucket i's list of points
	// each bucket_end[i] indexes one past the last element of bucket i's list of points
	bucket_begin = new device_vector<unsigned int>(prod(matrix_size_os+matrix_size_wrap));
	bucket_end = new device_vector<unsigned int>(prod(matrix_size_os+matrix_size_wrap));

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCUDA error in 'preprocess_generic_NFFT': %s. Quitting.", cudaGetErrorString(err)); fflush(stdout);
		exit(1);
	}

	// find the beginning of each bucket's list of points
	counting_iterator<unsigned int> search_begin(0);
	lower_bound(tuples_first.begin(), tuples_first.end(), search_begin,	search_begin + prod(matrix_size_os+matrix_size_wrap), bucket_begin->begin() );

	// find the end of each bucket's list of points
	upper_bound(tuples_first.begin(), tuples_first.end(), search_begin, search_begin + prod(matrix_size_os+matrix_size_wrap), bucket_end->begin() );
/*
	std::cout << std::endl << std::endl;
	for( int i=0; i<bucket_begin.size(); i++ )
	{
		unsigned int t = bucket_end[i];
		std::cout << "(" << t << "," << p.y << ")" << std::endl;
	}
*/
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'preprocess_generic_NFFT': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	successfully_preprocessed = true;
		
	return true;
}

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_iteration_plan< UINTd, FLOATd, TYPE >::preprocess_generic( unsigned int num_samples, FLOATd *traj )
{
	successfully_preprocessed = false;

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

	if(	traj_positions ){
		delete traj_positions;
		traj_positions = 0x0;
	}

	if(	tuples_last ){
		delete tuples_last;
		tuples_last = 0x0;
	}

	if(	bucket_begin ){
		delete bucket_begin;
		bucket_begin = 0x0;
	}

	if(	bucket_end ){
		delete bucket_end;
		bucket_end = 0x0;
	}

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'NFFT_iteration_plan::preprocess_generic': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Nullify radial only entries
	this->samples_per_projection = 0;
	this->projections_per_frame = 0;
	this->angular_offset = 0;
	this->frames_per_rotation = 0;
	this->interframe_rotation = 0;
	this->gc_factor = 0;
	this->total_projections_f = 0;
	this->number_of_strips_NFFT_H = 0;

	max_strips_per_domain_NFFT = 0;
	number_of_strips_NFFT = 0;
//	this->number_of_threads_NFFT_H = pre_NFFT_H->number_of_threads_NFFT_H;
	number_of_strips_NFFT_H = 0;
	stripsMapDevPtr_NFFT = 0x0;
	stripsDirDevPtr_NFFT = 0x0;
	stripOriginsDevPtr_NFFT = 0x0;
	stripLengthsDevPtr_NFFT = 0x0;
	domainsMapDevPtr_NFFT_H = 0x0;
	stripsMapDevPtr_NFFT_H = 0x0;
	stripsDevPtr_NFFT_H = 0x0;

	number_of_samples = num_samples;
	domain_count_samples = number_of_samples/domain_size_samples;

	NFFT_H_plan< UINTd, FLOATd, TYPE > *pre_NFFT_H = new NFFT_H_plan< UINTd, FLOATd, TYPE >( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_coils, W );
	pre_NFFT_H->preprocess_generic( num_samples, traj );

	if( !pre_NFFT_H->successfully_preprocessed ){
		NFFT_cleanup(&pre_NFFT_H);
		return false;
	}

	this->traj_positions = pre_NFFT_H->traj_positions;
	this->tuples_last = pre_NFFT_H->tuples_last;
	this->bucket_begin = pre_NFFT_H->bucket_begin;
	this->bucket_end = pre_NFFT_H->bucket_end;

	delete pre_NFFT_H;

	successfully_preprocessed = true;
		
	return true;
}

template <class FLOATd> __global__ void
noise_decorrelate_generic_kernel( unsigned int num_points, unsigned int num_coils, float shutter_radius_squared, cuFloatComplex *data, FLOATd *co, float *out_magnitudes_squared, float *out_count )
{
	const unsigned int sample_no = (blockIdx.x*blockDim.x+threadIdx.x);

	if( sample_no >= num_points )
		return;

	// Get sample position. 
	FLOATd sample_pos = co[sample_no];
	float dist_sq = dot(sample_pos, sample_pos);

	// Write to global memory

	if( dist_sq >= shutter_radius_squared ){
		for( unsigned int c=0; c<num_coils; c++ ){
			cuFloatComplex _data = data[c*num_points+sample_no];
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

template< class FLOATd > bool
noise_decorrelate_generic( unsigned int num_samples, unsigned int num_coils, float shutter_radius, cuFloatComplex *samples_DevPtr, FLOATd *trajectory_DevPtr )
{
	// First find noise magnitudes

	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)num_samples/(double)blockDim.x), 1, 1 );

	assert(num_samples>=blockDim.x);

	float *noise_modulus_sq, *point_count;

	cudaMalloc( (void**) &noise_modulus_sq, num_samples*num_coils*sizeof(float) );
	cudaMalloc( (void**) &point_count, num_samples*sizeof(float) );

	noise_decorrelate_generic_kernel<FLOATd><<< gridDim, blockDim >>>( num_samples, num_coils, shutter_radius*shutter_radius, samples_DevPtr, trajectory_DevPtr, noise_modulus_sq, point_count );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'noise_decorrelate_generic': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Then sum for each coil and normalize

	float *noise_variances = (float*) malloc(num_coils*sizeof(float));
	float used_points = cublasSasum( num_samples, point_count, 1 );
	if( used_points > 0.1f ){
		for( unsigned int c=0; c<num_coils; c++ ){
			noise_variances[c] = cublasSasum( num_samples, &noise_modulus_sq[c*num_samples], 1 );	
			noise_variances[c] /= used_points;
		}
	}
	else{
		printf("\nnoise_scaling: No points found within k-space shutter!\n");
		return false;
	}

	// And finally scale the samples values accordingly

	for( unsigned int c=0; c<num_coils; c++){
				//printf("\nNoise scaling with a factor of %f for coil %d", 1.0f/sqrtf(noise_variances[c]), c );
		cublasSscal( 2*num_samples, 1.0f/sqrtf(noise_variances[c]), (float*)(&samples_DevPtr[c*num_samples]), 1 );
	}

	free(noise_variances);
	cudaFree( noise_modulus_sq );
	cudaFree( point_count );

	return true;
}

// Set DCW
template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
bool set_dcw( PLAN<UINTd, FLOATd, TYPE> *plan, float *dcwDevPtr )
{
	cudaError_t err;

	if( plan->weights_DevPtr ){

		cudaFree( plan->weights_DevPtr );
		plan->weights_DevPtr = 0x0;

		err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nError freeing memory: 'set_dcw': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
			exit(1);
		}
	}

	cudaMalloc( (void**) &plan->weights_DevPtr, plan->number_of_samples*sizeof(float) );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nError allocating memory: 'set_dcw': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	cudaMemcpy( plan->weights_DevPtr, dcwDevPtr, plan->number_of_samples*sizeof(float), cudaMemcpyDeviceToDevice );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nError in memcpy: 'set_dcw': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	return true;
}

// Instantiation

template mr_recon::NFFT_H_plan<uint2, float2, 0>* preprocess_generic_NFFT( uint2, uint2, uint2, uint2, unsigned int, float, unsigned int, float2* );
template mr_recon::NFFT_H_plan<uint3, float3, 0>* preprocess_generic_NFFT( uint3, uint3, uint3, uint3, unsigned int, float, unsigned int, float3* );
template mr_recon::NFFT_H_plan<uint4, float4, 0>* preprocess_generic_NFFT( uint4, uint4, uint4, uint4, unsigned int, float, unsigned int, float4* );

template mr_recon::NFFT_iteration_plan<uint2, float2, 0>* preprocess_generic_NFFT( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, unsigned int, float2* );
template mr_recon::NFFT_iteration_plan<uint3, float3, 0>* preprocess_generic_NFFT( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, unsigned int, float3* );
template mr_recon::NFFT_iteration_plan<uint4, float4, 0>* preprocess_generic_NFFT( uint4, uint4, uint4, uint4, unsigned int, unsigned int, float, unsigned int, float4* );

template bool preprocess_generic_NFFT( mr_recon::NFFT_plan< uint2, float2, 0>*, unsigned int, float2* );
template bool preprocess_generic_NFFT( mr_recon::NFFT_plan< uint3, float3, 0>*, unsigned int, float3* );
template bool preprocess_generic_NFFT( mr_recon::NFFT_plan< uint4, float4, 0>*, unsigned int, float4* );

template bool preprocess_generic_NFFT( mr_recon::NFFT_H_plan< uint2, float2, 0>*, unsigned int, float2* );
template bool preprocess_generic_NFFT( mr_recon::NFFT_H_plan< uint3, float3, 0>*, unsigned int, float3* );
template bool preprocess_generic_NFFT( mr_recon::NFFT_H_plan< uint4, float4, 0>*, unsigned int, float4* );

template bool preprocess_generic_NFFT( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, unsigned int, float2* );
template bool preprocess_generic_NFFT( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, unsigned int, float3* );
template bool preprocess_generic_NFFT( mr_recon::NFFT_iteration_plan< uint4, float4, 0>*, unsigned int, float4* );

template bool noise_decorrelate_generic( unsigned int, unsigned int, float, cuFloatComplex*, float2* );
template bool noise_decorrelate_generic( unsigned int, unsigned int, float, cuFloatComplex*, float3* );
template bool noise_decorrelate_generic( unsigned int, unsigned int, float, cuFloatComplex*, float4* );

template bool set_dcw( mr_recon::NFFT_H_plan< uint2, float2, 0>*, float* );
template bool set_dcw( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, float* );
template bool set_dcw( mr_recon::NFFT_H_plan< uint3, float3, 0>*, float* );
template bool set_dcw( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, float* );