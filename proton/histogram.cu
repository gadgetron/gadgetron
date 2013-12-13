

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <thrust/inner_product.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>


#include "histogram.h"

texture<float, cudaTextureType1D, cudaReadModeElementType> texRef;

using namespace Gadgetron;


template<class T> __global__ void parzen_kernel(T* data,float * parzen,
		const float sigma, const float start,  const float end, const int bins, const int size){
	typedef typename realType<T>::type REAL;

	extern __shared__ T work_data[];

	int bin = threadIdx.x+blockIdx.y*blockDim.x;
	float res = 0;
	float val = bin*(start-end)/bins+start;
	float sigma2 = sigma*sigma;
	int idx = threadIdx.x+blockIdx.x*blockDim.x;

	if (idx < size)	work_data[threadIdx.x] = data[idx];
	__syncthreads();

	for (int i = 0; i < blockDim.x-(idx/size)*(idx%size); i++) {
		res += exp(-(val-work_data)*(val-work_data)/sigma2);
	}

	atomicAdd(&(parzen[bin]),res);
}


template<class T> boost::shared_ptr<cuNDArray<float> >  parzen_grid(cuNDArray<T> * data, int bins, float sigma, float& start, float & end ){

	int blockSize =128;
	dim3 gridSize((data->get_number_of_elements()-1)/blockSize+1,(bins-1)/blockSize+1);

	start = thrust::min_element(data->begin(),data->end());
	end = thrust::max_element(data->begin(),data->end());

	std::vector<size_t> pdims(1,bins);
	boost::shared_ptr<cuNDArray<float> > parzen(new cuNDArray<float>(&pdims));
	parzen_kernel<<<gridSize,blockSize,blockSize*sizeof(T)>>>(data->get_data_ptr(),parzen->get_data_ptr(),sigma,start,end,bins,data->get_number_of_elements());

	return parzen;
}


template<class T> __global__ void entropy_kernel(T* data,float * parzen,float* result, int size, float start, float stop){
	typedef typename realType<T>::type REAL;

	extern __shared__ float work_data[];
	int idx = threadIdx.x+blockIdx.x*blockDim.x;
	if (idx < size){
		 float p = tex1D(texRef,(float(data[idx])-start)/(stop-start));
		 work_data[threadIdx.x] = p*log2(p);
		__syncthreads();
		for(int offset = blockDim.x / 2; offset > 0; offset /= 2)
		  {
		    if(threadIdx.x < offset)
		    {
		     work_data[threadIdx.x] += work_data[threadIdx.x + offset];
		    }
		    __syncthreads();
		  }


	}
	if (threadIdx.x == 0) atomicAdd(result,work_data[0]);


}


template<class T> float entropy(cuNDArray<T> * data, int bins, float sigma){
	float start,end;
	boost::shared_ptr<cuNDArray<float> > grid = parzen_grid(data,bins,sigma,start, end);
	int blockSize =128;
	int gridSize = (data->get_number_of_elements()-1)/blockSize+1;

	texRef.filterMode = cudaFilterModeLinear;
	texRef.normalized = true;
	size_t offset;
	cudaBindTexture(&offset,texRef,grid->get_data_ptr(),grid->get_number_of_elements());
	thrust::device_vector<float> result_vec(1,0.0f);
	float* dev_result = thrust::raw_pointer_cast( &result_vec[0] );
	entropy_kernel<<<gridSize,blockSize,blockSize*sizeof(float)>>>(data->get_data_ptr(),grid->get_data_ptr(),sigma,start,end,bins,data->get_number_of_elements());

	return result_vec[0];
}

