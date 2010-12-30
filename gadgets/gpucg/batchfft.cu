/**
 * @file
 * Implementation of batched 2D FFTs with CUDA
 *
 * @author Jim Hardwick
 *
 * @note
 * This code was posted to nVidia's public CUDA forum as thanks to the CUDA
 * forum membership for their assistance. No copyright was asserted and no
 * license terms were imposed in the post. The post can be found at
 * http://forums.nvidia.com/index.php?showtopic=34241
 */

/**
 * @note
 * Modified and extended to suit the mr_recon project...
 */

#include <cufft.h>

#include "batchfft.hcu"
#include "uint_util.hcu"

#define BLOCK_DIM 16 ///< Thread block dimension. (TSS: I think think should be 'half_warp_size' "instead"?)

texture<float2, 1> permute_tex;

// Choose which "3D" transpose
void prepare_xy_transpose( batchfftHandle *plan );
//void prepare_yz_transpose( batchfftHandle *plan );
void prepare_xz_transpose( batchfftHandle *plan );

/**
 * Fast matrix transpose kernel
 *
 * Uses shared memory to coalesce global memory reads and writes, improving performance.
 *
 * @param idata Input batch of matrices
 * @param odata Output buffer for transposed batch of matrices, must be different than idata
 * @param width Width of each matrix, must be a multiple of BLOCK_DIM
 * @param height Height of each matrix, must be a multiple of BLOCK_DIM
 * @param num Matrices in the batch
 */
__global__ void transpose_xy(float2* idata, float2* odata, int width, int height, int num);
__global__ void transpose_xz(float2* idata, float2* odata, int width, int height, int num_y, int num_batch);
__global__ void transpose_zx(float2* idata, float2* odata, int width, int height, int num_y, int num_batch);

// Permute kernel ("3D transpose")
__global__ void FFT_permute_kernel( cuFloatComplex* data_in, cuFloatComplex* data_out, uint4 dim, unsigned int num_shifts );

////////////////////////////////////////////////////////////////////////////////
cufftResult batchfftPlan2d(batchfftHandle* plan, int nx, int ny, cufftType type, int batch)
{
	if(type != CUFFT_C2C)
		return CUFFT_INVALID_TYPE;

	if((nx % BLOCK_DIM) != 0)
		return CUFFT_INVALID_SIZE;

	if((ny % BLOCK_DIM) != 0)
		return CUFFT_INVALID_SIZE;

	// Swap nx and ny so they correspoind to the 2D CUFFT API.
	plan->nx = ny;
	plan->ny = nx;
	plan->nz = 1;
	plan->type = type;
	plan->batch = batch;

	plan->transpose_threads.x = BLOCK_DIM;
	plan->transpose_threads.y = BLOCK_DIM;
	plan->transpose_threads.z = 1;
	plan->transpose_grid.x = plan->nx / BLOCK_DIM;
	plan->transpose_grid.y = plan->ny / BLOCK_DIM;
	plan->transpose_grid.z = 1;
	plan->transpose_back_grid.x = plan->ny / BLOCK_DIM;
	plan->transpose_back_grid.y = plan->nx / BLOCK_DIM;
	plan->transpose_back_grid.z = 1;

	cufftResult ret = CUFFT_SUCCESS;
	cudaError_t cudaret = cudaSuccess;
	
	cudaret = cudaMalloc(&(plan->temp), plan->nx * plan->ny * plan->batch * sizeof(float2));
	if(cudaret != cudaSuccess)
		return CUFFT_ALLOC_FAILED;

	ret = cufftPlan1d(&(plan->rowplan), plan->nx, plan->type, plan->ny * plan->batch);
	if(ret != CUFFT_SUCCESS)
	{
		cudaret = cudaFree(plan->temp);
		if(cudaret != cudaSuccess){
			printf("\n'cudaFree' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}
		plan->temp = NULL;
		return ret;
	}

	ret = cufftPlan1d(&(plan->colplan), plan->ny, plan->type, plan->nx * plan->batch);
	if(ret != CUFFT_SUCCESS)
	{
		cudaret = cudaFree(plan->temp);
		if(cudaret != cudaSuccess){
			printf("\n'cudaFree' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
		}
		plan->temp = NULL;
		
		printf("\n 'cufftPlan1d' failed!."); fflush(stdout);
		
		cufftResult res2 = cufftDestroy(plan->rowplan);
		if(res2 != CUFFT_SUCCESS)
		{
			printf("\n'cufftDestroy' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		return ret;
	}

	return CUFFT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
cufftResult batchfftPlan3d(batchfftHandle* plan, int nx, int ny, int nz, cufftType type, int batch)
{
	if(type != CUFFT_C2C)
		return CUFFT_INVALID_TYPE;

	if((nx % BLOCK_DIM) != 0){
		printf("\nbatchFFT: all dimensions (current: %d) must be a multiple of the BLOCK_DIM (%d)", nx, BLOCK_DIM );
		return CUFFT_INVALID_SIZE;
	}

	if((ny % BLOCK_DIM) != 0){
		printf("\nbatchFFT: all dimensions (current: %d) must be a multiple of the BLOCK_DIM (%d)", ny, BLOCK_DIM );
		return CUFFT_INVALID_SIZE;
	}

	if((nz % BLOCK_DIM) != 0){
		printf("\nbatchFFT: all dimensions (current: %d) must be a multiple of the BLOCK_DIM (%d)", nz, BLOCK_DIM );
		return CUFFT_INVALID_SIZE;
	}

	// Swap nx, ny, and nz so they correspoind to the 3D CUFFT API.
	plan->nx = nz;
	plan->ny = ny;
	plan->nz = nx;
	plan->type = type;
	plan->batch = batch;

	plan->transpose_threads.x = BLOCK_DIM;
	plan->transpose_threads.y = BLOCK_DIM;
	plan->transpose_threads.z = 1;

	cufftResult ret = CUFFT_SUCCESS;
	cudaError_t cudaret = cudaSuccess;
	
	cudaret = cudaMalloc(&(plan->temp), plan->nx * plan->ny * plan->nz * plan->batch * sizeof(float2));
	if(cudaret != cudaSuccess)
		return CUFFT_ALLOC_FAILED;

	ret = cufftPlan1d(&(plan->rowplan), plan->nx, plan->type, plan->ny * plan->nz * plan->batch);
	if(ret != CUFFT_SUCCESS)
	{
		cudaret = cudaFree(plan->temp);
		if(cudaret != cudaSuccess){
			printf("\n'cudaFree' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		plan->temp = NULL;
		return ret;
	}

	ret = cufftPlan1d(&(plan->colplan), plan->ny, plan->type, plan->nx * plan->nz * plan->batch);
	if(ret != CUFFT_SUCCESS)
	{
		cudaret = cudaFree(plan->temp);
		if(cudaret != cudaSuccess){
			printf("\n'cudaFree' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		printf("\n 'cufftPlan1d' failed!."); fflush(stdout);

		plan->temp = NULL;
		cufftResult res2 = cufftDestroy(plan->rowplan);

		if(res2 != CUFFT_SUCCESS)
		{
			printf("\n'cufftDestroy' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		return ret;
	}

	ret = cufftPlan1d(&(plan->layerplan), plan->nz, plan->type, plan->nx * plan->ny * plan->batch);
	if(ret != CUFFT_SUCCESS)
	{
		cudaret = cudaFree(plan->temp);
		if(cudaret != cudaSuccess){
			printf("\n'cudaFree' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		plan->temp = NULL;

		cufftResult res2 = cufftDestroy(plan->rowplan);
		if(res2 != CUFFT_SUCCESS)
		{
			printf("\n'cufftDestroy' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		res2 = cufftDestroy(plan->colplan);
		if(res2 != CUFFT_SUCCESS)
		{
			printf("\n'cufftDestroy' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		return ret;
	}

	return CUFFT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
cufftResult batchfftDestroy(batchfftHandle* plan)
{
	if(plan->temp != NULL)
	{
		cufftResult res = cufftDestroy(plan->rowplan);
		if(res != CUFFT_SUCCESS)
		{
			printf("\n'cufftDestroy' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		res = cufftDestroy(plan->colplan);
		if(res != CUFFT_SUCCESS)
		{
			printf("\n'cufftDestroy' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		if( plan->nz != 1 ){
			res = cufftDestroy(plan->layerplan);
			if(res != CUFFT_SUCCESS)
			{
			printf("\n'cufftDestroy' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
				exit(1);
			}
		}

		cudaError_t cudaret = cudaFree(plan->temp);
		if(cudaret != cudaSuccess){
			printf("\n'cudaFree' failed: %s. Quitting.", cudaGetErrorString(cudaGetLastError())); fflush(stdout);
			exit(1);
		}

		plan->temp = NULL;
	}

	return CUFFT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
cufftResult batchfftExecute(batchfftHandle plan, void* idata, void* odata, int sign)
{
	if( plan.nz == 1 )
	{
		// 2D FFT

		cufftResult cufftret = CUFFT_SUCCESS;
		cudaError_t cudaret = cudaSuccess;

		// Transform rows
		cufftret = cufftExecC2C(plan.rowplan, (cuFloatComplex*)idata, (cuFloatComplex*)odata, sign);
		if(cufftret != CUFFT_SUCCESS)
			return cufftret;

		// Transpose
		transpose_xy<<< plan.transpose_grid, plan.transpose_threads >>>((float2*)odata, (float2*)plan.temp, plan.nx, plan.ny, plan.batch);
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;

		// Transform columns
		cufftret = cufftExecC2C(plan.colplan, (cuFloatComplex*)plan.temp, (cuFloatComplex*)plan.temp, sign);
		if(cufftret != CUFFT_SUCCESS)
			return cufftret;

		// Transpose back
		transpose_xy<<< plan.transpose_back_grid, plan.transpose_threads >>>((float2*)plan.temp, (float2*)odata, plan.ny, plan.nx, plan.batch);
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;
	}
	else{

		// 3D FFT

		cufftResult cufftret = CUFFT_SUCCESS;
		cudaError_t cudaret = cudaSuccess;

		dim3 blockDim(512,1, 1);
		dim3 gridDim((int)ceil((double)(plan.nx*plan.ny*plan.nz*plan.batch)/512.0),1,1);

		const uint4 dims = make_uint4( plan.nx, plan.ny, plan.nz, plan.batch );
		const cudaChannelFormatDesc desc_float2 = cudaCreateChannelDesc<float2>();

		// Transform rows
		cufftret = cufftExecC2C(plan.rowplan, (cuFloatComplex*)idata, (cuFloatComplex*)odata, sign);
		if(cufftret != CUFFT_SUCCESS)
			return cufftret;

//		fft_permute_radix4( plan.nx, dims, 0, plan.ny*plan.nz*plan.batch, idata, odata, sign );
/*
		// Permute once "to the left"
		cudaBindTexture( 0, &permute_tex, (float2*) odata, &desc_float2, prod(dims)*sizeof(float2));
		FFT_permute_kernel<<< gridDim, blockDim >>>( (cuFloatComplex*) odata, (cuFloatComplex*) plan.temp, make_uint4(plan.nx, plan.ny, plan.nz, plan.batch), 1 );
		cudaUnbindTexture( permute_tex );
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;
*/

		// Transpose
		prepare_xy_transpose( &plan );
		transpose_xy<<< plan.transpose_grid, plan.transpose_threads >>>((float2*)odata, (float2*)plan.temp, plan.nx, plan.ny, plan.nz*plan.batch);
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;

		// Transform columns
		cufftret = cufftExecC2C(plan.colplan, (cuFloatComplex*)plan.temp, (cuFloatComplex*) plan.temp, sign);
		if(cufftret != CUFFT_SUCCESS)
			return cufftret;

		// Transpose back
		prepare_xy_transpose( &plan );
		transpose_xy<<< plan.transpose_back_grid, plan.transpose_threads >>>((float2*)plan.temp, (float2*)odata, plan.ny, plan.nx, plan.nz*plan.batch);
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;

//		fft_permute_radix4( plan.ny, dims, 1, plan.nx*plan.nz*plan.batch, odata, odata, sign );
/*
		// Permute once "to the left"
		cudaBindTexture( 0, &permute_tex, (float2*) plan.temp, &desc_float2, prod(dims)*sizeof(float2));
		FFT_permute_kernel<<< gridDim, blockDim >>>( (cuFloatComplex*) plan.temp, (cuFloatComplex*) odata, make_uint4(plan.ny, plan.nz, plan.batch, plan.nx), 1 );
		cudaUnbindTexture( permute_tex );
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;
*/

		// Transpose
		prepare_xz_transpose( &plan );
		transpose_xz<<< plan.transpose_grid, plan.transpose_threads >>>((float2*)odata, (float2*)plan.temp, plan.nx, plan.nz, plan.ny, plan.batch);
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;

		// Transform layers
		cufftret = cufftExecC2C(plan.layerplan, (cuFloatComplex*) plan.temp, (cuFloatComplex*)plan.temp, sign);
		if(cufftret != CUFFT_SUCCESS)
			return cufftret;
		
		// Transpose back
		prepare_xz_transpose( &plan );
		transpose_zx<<< plan.transpose_back_grid, plan.transpose_threads >>>((float2*)plan.temp, (float2*)odata, plan.nz, plan.nx, plan.ny, plan.batch);
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;

//		fft_permute_radix2( plan.nz, dims, 2, plan.nx*plan.ny*plan.batch, odata, odata, sign );
/*
		// Permute one final time to complete the loop
		cudaBindTexture( 0, &permute_tex, (float2*) plan.temp, &desc_float2, prod(dims)*sizeof(float2));
		FFT_permute_kernel<<< gridDim, blockDim >>>( (cuFloatComplex*) plan.temp, (cuFloatComplex*) odata, make_uint4(plan.nz, plan.batch, plan.nx, plan.ny), 2 );
		cudaUnbindTexture( permute_tex );
		cudaret = cudaGetLastError();
		if(cudaret != cudaSuccess)
			return CUFFT_EXEC_FAILED;
*/
	}

	return CUFFT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
__global__ void transpose_xy(float2* idata, float2* odata, int width, int height, int num)
{
	// To prevent shared memory bank confilcts:
	// - Load each component into a different array. Since the array size is a
	//   multiple of the number of banks (16), each thread reads x and y from
	//   the same bank. If a single float2 array is used, thread n would read
	//   x and y from banks n and n+1, and thread n+8 would read values from the
	//   same banks - causing a bank conflict.
	// - Use BLOCK_DIM+1 as the x size of the array. This way each row of the
	//   array starts in a different bank - so reading from shared memory
	//   doesn't cause bank conflicts when writing the transpose out to global
	//   memory.
	__shared__ float blockx[(BLOCK_DIM+1)*BLOCK_DIM];
	__shared__ float blocky[(BLOCK_DIM+1)*BLOCK_DIM];

	unsigned int xBlock = BLOCK_DIM * blockIdx.x;
	unsigned int yBlock = BLOCK_DIM * blockIdx.y;
	unsigned int xIndex = xBlock + threadIdx.x;
	unsigned int yIndex = yBlock + threadIdx.y;
	unsigned int size = width * height;

	for(int n = 0; n < num; ++n){
		unsigned int index_in = n * size + width * yIndex + xIndex;
		unsigned int index_block = threadIdx.y * (BLOCK_DIM+1) + threadIdx.x;

		float2 in = idata[index_in];
		blockx[index_block] = in.x;
		blocky[index_block] = in.y;

		unsigned int index_transpose = threadIdx.x * (BLOCK_DIM+1) + threadIdx.y;
		unsigned int index_out = n * size + height * (xBlock + threadIdx.y) + yBlock + threadIdx.x;

		__syncthreads();

		float2 out = make_float2( blockx[index_transpose], blocky[index_transpose] );
		odata[index_out] = out;
	}
} 

__global__ void transpose_xz(float2* idata, float2* odata, int width, int height, int num_y, int num_batch)
{
	// As above, but the 'y' direction elements are now a slice rather than a line apart

	__shared__ float blockx[(BLOCK_DIM+1)*BLOCK_DIM];
	__shared__ float blocky[(BLOCK_DIM+1)*BLOCK_DIM];

	unsigned int xBlock = BLOCK_DIM * blockIdx.x;
	unsigned int yBlock = BLOCK_DIM * blockIdx.y;
	unsigned int xIndex = xBlock + threadIdx.x;
	unsigned int yIndex = yBlock + threadIdx.y;
	unsigned int size = width * height;
	unsigned int orig_size = width*num_y;
	unsigned int volume = size*num_y;

	for( int n_b = 0; n_b < num_batch; ++n_b ){
		for( int n_y = 0; n_y < num_y; ++n_y ){
//		unsigned int index_in = n * size + width * yIndex + xIndex;
			unsigned int index_in = n_b * volume + n_y * width + orig_size * yIndex + xIndex;
			unsigned int index_block = threadIdx.y * (BLOCK_DIM+1) + threadIdx.x;

			float2 in = idata[index_in];
			blockx[index_block] = in.x;
			blocky[index_block] = in.y;

			unsigned int index_transpose = threadIdx.x * (BLOCK_DIM+1) + threadIdx.y;
//			unsigned int index_out = n * size + height * (xBlock + threadIdx.y) + yBlock + threadIdx.x;
			unsigned int index_out = n_b * volume + n_y * size + height * (xBlock + threadIdx.y) + yBlock + threadIdx.x;

			__syncthreads();

			float2 out = make_float2( blockx[index_transpose], blocky[index_transpose] );
			odata[index_out] = out;
		}
	}
}

__global__ void transpose_zx(float2* idata, float2* odata, int width, int height, int num_y, int num_batch)
{
	// As above, but the 'y' direction elements are now a slice rather than a line apart

	__shared__ float blockx[(BLOCK_DIM+1)*BLOCK_DIM];
	__shared__ float blocky[(BLOCK_DIM+1)*BLOCK_DIM];

	unsigned int xBlock = BLOCK_DIM * blockIdx.x;
	unsigned int yBlock = BLOCK_DIM * blockIdx.y;
	unsigned int xIndex = xBlock + threadIdx.x;
	unsigned int yIndex = yBlock + threadIdx.y;
	unsigned int size = width * height;
	unsigned int orig_size = height*num_y;
	unsigned int volume = size*num_y;

	for( int n_b = 0; n_b < num_batch; ++n_b ){
		for( int n_y = 0; n_y < num_y; ++n_y ){
//		unsigned int index_in = n * size + width * yIndex + xIndex;
			unsigned int index_in = n_b * volume + n_y * size + width * yIndex + xIndex;
			unsigned int index_block = threadIdx.y * (BLOCK_DIM+1) + threadIdx.x;

			float2 in = idata[index_in];
			blockx[index_block] = in.x;
			blocky[index_block] = in.y;

			unsigned int index_transpose = threadIdx.x * (BLOCK_DIM+1) + threadIdx.y;
//			unsigned int index_out = n * size + height * (xBlock + threadIdx.y) + yBlock + threadIdx.x;
			unsigned int index_out = n_b * volume + n_y * height + orig_size * xBlock + yBlock + orig_size * threadIdx.y + threadIdx.x;

			__syncthreads();

			float2 out = make_float2( blockx[index_transpose], blocky[index_transpose] );
			odata[index_out] = out;
		}
	}
}

/*
__global__ void FFT_permute_kernel( cuFloatComplex* data_in, cuFloatComplex* data_out, uint4 dim, unsigned int num_shifts )
{	
	//This is the current pixel number
	unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	// Do nothing if index is out of range
	if( idx < prod(dim) ){

		uint4 new_dim;
		uint4 co = idx_to_co(idx, dim);

		co = shift_down( co, num_shifts );
		new_dim = shift_down( dim, num_shifts );

		data_out[co_to_idx(co,new_dim)] = data_in[idx];
	}
}
*/


__global__ void FFT_permute_kernel( cuFloatComplex* data_in, cuFloatComplex* data_out, uint4 dim, unsigned int num_shifts )
{	
	//This is the current pixel number
	unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	// Do nothing if index is out of range
	if( idx < prod(dim) ){

		uint4 new_dim = shift_down( dim, num_shifts );
		uint4 co = idx_to_co(idx, new_dim);
		co = shift_up( co, num_shifts );
		
//		data_out[idx] = data_in[co_to_idx(co,dim)];
		data_out[idx] = tex1Dfetch( permute_tex, co_to_idx(co,dim) );
	}
}


void prepare_xy_transpose( batchfftHandle *plan )
{
	plan->transpose_grid.x = plan->nx / BLOCK_DIM;
	plan->transpose_grid.y = plan->ny / BLOCK_DIM;
	plan->transpose_grid.z = 1;
	plan->transpose_back_grid.x = plan->ny / BLOCK_DIM;
	plan->transpose_back_grid.y = plan->nx / BLOCK_DIM;
	plan->transpose_back_grid.z = 1;
}
/*
void prepare_yz_transpose( batchfftHandle *plan )
{
	plan->transpose_grid.x = plan->ny / BLOCK_DIM;
	plan->transpose_grid.y = plan->nz / BLOCK_DIM;
	plan->transpose_grid.z = 1;
	plan->transpose_back_grid.x = plan->nz / BLOCK_DIM;
	plan->transpose_back_grid.y = plan->ny / BLOCK_DIM;
	plan->transpose_back_grid.z = 1;
}
*/
void prepare_xz_transpose( batchfftHandle *plan )
{
	plan->transpose_grid.x = plan->nx / BLOCK_DIM;
	plan->transpose_grid.y = plan->nz / BLOCK_DIM;
	plan->transpose_grid.z = 1;
	plan->transpose_back_grid.x = plan->nz / BLOCK_DIM;
	plan->transpose_back_grid.y = plan->nx / BLOCK_DIM;
	plan->transpose_back_grid.z = 1;
}
