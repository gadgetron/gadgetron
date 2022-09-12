#pragma once
#include "vector_td.h"
#include "vector_td_utilities.h"
#include "setup_grid.h"
#include "cuSparseMatrix.h"
#include <thrust/fill.h>
#include "ConvolverNC2C_standard.cuh"
namespace Gadgetron {

template<unsigned int N> struct iteration_counter{};

template<class T, unsigned int D, unsigned int N, template<class, unsigned int> class K>
__device__ __inline__
realType_t<T> ndim_loop(const vector_td<realType_t<T>,D> point, unsigned int & loop_counter,
	T * __restrict__ weights,
		int * __restrict__ column_indices,
		vector_td<int, D> & grid_point, const vector_td<int,D> & image_dims,
		const ConvolutionKernel<realType_t<T>, D, K>* kernel,
		iteration_counter<N>)
{
	realType_t<T> wsum = 0;
	realType_t<T> radius = kernel->get_radius();
	for (int i = ::ceil(point[N] - radius); i <= ::floor(point[N] + radius); i++)
	{
		grid_point[N] = i;
		wsum += ndim_loop(point,loop_counter,weights,column_indices,grid_point,
			image_dims,kernel,iteration_counter<N-1>());
	}
	return wsum;
}

template<class T, unsigned int D, template<class, unsigned int> class K>
__device__ __inline__
realType_t<T> ndim_loop(const vector_td<realType_t<T>,D> point, unsigned int & loop_counter,
		T * __restrict__ weights,
		int * __restrict__ column_indices,
		vector_td<int, D> & grid_point, const vector_td<int,D> & image_dims,
		const ConvolutionKernel<realType_t<T>, D, K>* kernel,
		iteration_counter<0>)
{
	realType_t<T> wsum =0;
	realType_t<T> radius = kernel->get_radius();
	for (int i = ::ceil(point[0] - radius); i <= ::floor(point[0] + radius); i++)
	{
		grid_point[0] = i;
		realType_t<T> weight = kernel->get(abs(point-vector_td<realType_t<T>,D>(grid_point)));
		weights[loop_counter] = weight;
		//column_indices[loop_counter] = co_to_idx(grid_point%image_dims,image_dims);
		column_indices[loop_counter] = co_to_idx((grid_point+image_dims)%image_dims,image_dims);
		loop_counter++;
		wsum += weight;
	}
	return wsum;
}

template<class T> __device__ void index_sort(T* values, int* indices, int nvals){

	//Insertion sort, as we anticipate stuff to be mostly sorted. Might be faster to just sort the entire thing in a separate kernel
	for (int i = 0; i < nvals; i++)
	{
		int index = indices[i];
		T val = values[i];
		int j = i;
		while (j > 0 && indices[j-1] > index){
			indices[j] = indices[j-1];
			values[j] = values[j-1];
			j--;
		}
		indices[j] = index;
		values[j] = val;
	}
}

/**
 *
 * @param points non-cartesian points on which to do the NFFT. Size tot_size
 * @param offsets Array containing the offsets at which the rows should be stored. Size tot_size+1
 * @param weights Output array which will contain the values of the sparse matrix.
 * @param column_indices Output array containing the column indices of the sparse matrix
 * @param tot_size Total number of points
 */
template<class T, unsigned int D, template<class, unsigned int> class K>
__global__
void make_conv_matrix_kernel(const vector_td<realType_t<T>,D> * __restrict__ points,
	const int * __restrict__ offsets, T * __restrict__ weights,
	int * __restrict__ column_indices, const vector_td<int,D> image_dims,
	unsigned int tot_size,
	const ConvolutionKernel<realType_t<T>, D, K>* kernel
	)
{
	const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < tot_size){
		vector_td<realType_t<T>,D> p = points[idx];
		const int offset = offsets[idx];


		T * local_weight = weights+offset;
		int * local_column_indices = column_indices+offset;
		unsigned int loop_counter = 0;
		vector_td<int,D> grid_point;
		realType_t<T> wsum = ndim_loop(p,loop_counter, local_weight,local_column_indices,
			grid_point,image_dims,kernel,iteration_counter<D-1>());
		realType_t<T> inv_wsum = 1.0/wsum;
		for (unsigned int i = offset; i < offsets[idx+1]; i++){
			weights[i] *= inv_wsum;
		}
		index_sort(local_weight,local_column_indices,offsets[idx+1]-offset);

	}
}

template<class T>
void check_csrMatrix(cuCsrMatrix<T > &matrix)
{

	if (matrix.csrRow.size() != matrix.m+1){
		throw std::runtime_error("Malformed CSR matrix: length of CSR vector does not match matrix size");
	}

	if (matrix.csrColdnd.size() != matrix.nnz){
		throw std::runtime_error("Malformed CSR matrix: length of column indices vector does not match number of non-zero elements");
	}
	if (matrix.data.size() != matrix.nnz ){
		throw std::runtime_error("Malformed CSR matrix: length of data vector does not match number of non-zero elements");
	}

	int min_ind = *thrust::min_element(matrix.csrColdnd.begin(),matrix.csrColdnd.end());
	int max_ind = *thrust::max_element(matrix.csrColdnd.begin(),matrix.csrColdnd.end());

	if (min_ind < 0 || max_ind > matrix.n){
		std::stringstream ss;
		ss << "Malformed CSR matrix: column indices vector contains illegal values. Min " << min_ind<< " max " << max_ind;
		throw std::runtime_error(ss.str());
	}
	int min_row = *thrust::min_element(matrix.csrRow.begin(),matrix.csrRow.end());
	int max_row = *thrust::max_element(matrix.csrRow.begin(),matrix.csrRow.end());

	if (min_row < 0 || max_row != matrix.nnz){
		throw std::runtime_error("Malformed CSR matrix: CSR vector conains illegal values");
	}


	if (isnan(abs(thrust::reduce(matrix.data.begin(),matrix.data.end()))))
		throw std::runtime_error("Matrix contains NaN");


}

template<class T, unsigned int D, template<class, unsigned int> class K>
cuCsrMatrix<T> make_conv_matrix(
	const thrust::device_vector<vector_td<realType_t<T>,D>> &points, 
	const vector_td<size_t,D>& image_dims, 
	const ConvolutionKernel<realType_t<T>, D, K>* kernel)
{
	auto csrRow = thrust::device_vector<int>(points.size()+1);
	csrRow[0] = 0;
	CHECK_FOR_CUDA_ERROR();

	realType_t<T> radius = kernel->get_radius();
	{
		thrust::device_vector<int> c_p_s(points.size());
		thrust::transform(points.begin(), points.end(), c_p_s.begin(),
			compute_num_cells_per_sample<realType_t<T>,D>(kernel->get_radius()));

		thrust::inclusive_scan( c_p_s.begin(), c_p_s.end(), csrRow.begin()+1, thrust::plus<int>()); // prefix sum


	}
	unsigned int num_pairs = csrRow.back();
	//cuNDArray<int> row_indices(ind_dims);
	auto csrColdnd = thrust::device_vector<int>(num_pairs);
	auto data = thrust::device_vector<T >(num_pairs);
	//cuNDArray<T > values(ind_dims);


	dim3 dimBlock;
	dim3 dimGrid;
	setup_grid(points.size(),&dimBlock,&dimGrid);

	make_conv_matrix_kernel<<<dimGrid,dimBlock>>>(
		thrust::raw_pointer_cast(points.data()),
		thrust::raw_pointer_cast(csrRow.data()),
		thrust::raw_pointer_cast(data.data()),
		thrust::raw_pointer_cast(csrColdnd.data()),
		vector_td<int,D>(image_dims), points.size(),kernel);
	cudaDeviceSynchronize();
	CHECK_FOR_CUDA_ERROR();

	cuCsrMatrix<T> matrix(prod(image_dims), points.size(),std::move(csrRow),std::move(csrColdnd),std::move(data));
	return matrix;
}

}


