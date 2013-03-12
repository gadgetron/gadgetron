#include "cuVector_td_test_kernels.h"
#include "check_CUDA.h"
#include "vector_td_utilities.h"
#include "cuNDArray.h"
#include "cudaDeviceManager.h"
#include "thrust/device_vector.h"


using namespace Gadgetron;
template<class T, unsigned int D> __global__ void abs_kernel(vector_td<T,D>* data, unsigned int size){
	 const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	 if (idx < size) data[idx] = abs(data[idx]);
}


template<class T, unsigned int D> void Gadgetron::test_abs(cuNDArray< vector_td<T,D> >* data){

	dim3 dimBlock(std::min(cudaDeviceManager::Instance()->max_griddim(),(int)data->get_number_of_elements()));
	dim3 dimGrid((dimBlock.x-1)/data->get_number_of_elements()+1);
	abs_kernel<<<dimGrid,dimBlock>>>(data->get_data_ptr(),data->get_number_of_elements());
	cudaThreadSynchronize();
	CHECK_FOR_CUDA_ERROR();
}


template<typename T, unsigned D>
struct test_norm_functor : public thrust::unary_function<T,vector_td<T,D> >
{
 __host__ __device__ T operator()(const vector_td<T,D> &x) const {return norm(x);}
};
template<class T, unsigned int D> thrust::device_vector<T> Gadgetron::test_norm(cuNDArray< vector_td<T,D> >* data){

	thrust::device_vector<T> out(data->get_number_of_elements());
	thrust::transform(data->begin(),data->end(),out.begin(),test_norm_functor<T,D>());
	cudaThreadSynchronize();
	CHECK_FOR_CUDA_ERROR();
	return out;
}



template<typename T, unsigned D>
struct test_min_functor : public thrust::unary_function<T,vector_td<T,D> >
{
 __host__ __device__ T operator()(const vector_td<T,D> &x) const {return min(x);}
};
template<class T, unsigned int D> thrust::device_vector<T> Gadgetron::test_min(cuNDArray< vector_td<T,D> >* data){

	thrust::device_vector<T> out(data->get_number_of_elements());
	thrust::transform(data->begin(),data->end(),out.begin(),test_min_functor<T,D>());
	cudaThreadSynchronize();
	CHECK_FOR_CUDA_ERROR();
	return out;
}


template<typename T, unsigned D>
struct test_max_functor : public thrust::unary_function<T,vector_td<T,D> >
{
 __host__ __device__ T operator()(const vector_td<T,D> &x) const {return max(x);}
};
template<class T, unsigned int D> thrust::device_vector<T> Gadgetron::test_max(cuNDArray< vector_td<T,D> >* data){

	thrust::device_vector<T> out(data->get_number_of_elements());
	thrust::transform(data->begin(),data->end(),out.begin(),test_max_functor<T,D>());
	cudaThreadSynchronize();
	CHECK_FOR_CUDA_ERROR();
	return out;
}


template void Gadgetron::test_abs<float,1>(cuNDArray< vector_td<float,1> > *);
template void Gadgetron::test_abs<float,2>(cuNDArray< vector_td<float,2> > *);
template  void Gadgetron::test_abs<float,3>(cuNDArray< vector_td<float,3> > *);
template  void Gadgetron::test_abs<float,4>(cuNDArray< vector_td<float,4> > *);

template  void Gadgetron::test_abs<double,1>(cuNDArray< vector_td<double,1> > *);
template void Gadgetron::test_abs<double,2>(cuNDArray< vector_td<double,2> > *);
template void Gadgetron::test_abs<double,3>(cuNDArray< vector_td<double,3> > *);
template void Gadgetron::test_abs<double,4>(cuNDArray< vector_td<double,4> > *);


template thrust::device_vector<float> Gadgetron::test_norm<float,1>(cuNDArray< vector_td<float,1> > *);
template thrust::device_vector<float> Gadgetron::test_norm<float,2>(cuNDArray< vector_td<float,2> > *);
template  thrust::device_vector<float> Gadgetron::test_norm<float,3>(cuNDArray< vector_td<float,3> > *);
template  thrust::device_vector<float> Gadgetron::test_norm<float,4>(cuNDArray< vector_td<float,4> > *);

template  thrust::device_vector<double> Gadgetron::test_norm<double,1>(cuNDArray< vector_td<double,1> > *);
template thrust::device_vector<double> Gadgetron::test_norm<double,2>(cuNDArray< vector_td<double,2> > *);
template thrust::device_vector<double> Gadgetron::test_norm<double,3>(cuNDArray< vector_td<double,3> > *);
template thrust::device_vector<double> Gadgetron::test_norm<double,4>(cuNDArray< vector_td<double,4> > *);


template thrust::device_vector<float> Gadgetron::test_min<float,1>(cuNDArray< vector_td<float,1> > *);
template thrust::device_vector<float> Gadgetron::test_min<float,2>(cuNDArray< vector_td<float,2> > *);
template  thrust::device_vector<float> Gadgetron::test_min<float,3>(cuNDArray< vector_td<float,3> > *);
template  thrust::device_vector<float> Gadgetron::test_min<float,4>(cuNDArray< vector_td<float,4> > *);

template  thrust::device_vector<double> Gadgetron::test_min<double,1>(cuNDArray< vector_td<double,1> > *);
template thrust::device_vector<double> Gadgetron::test_min<double,2>(cuNDArray< vector_td<double,2> > *);
template thrust::device_vector<double> Gadgetron::test_min<double,3>(cuNDArray< vector_td<double,3> > *);
template thrust::device_vector<double> Gadgetron::test_min<double,4>(cuNDArray< vector_td<double,4> > *);


template thrust::device_vector<float> Gadgetron::test_max<float,1>(cuNDArray< vector_td<float,1> > *);
template thrust::device_vector<float> Gadgetron::test_max<float,2>(cuNDArray< vector_td<float,2> > *);
template  thrust::device_vector<float> Gadgetron::test_max<float,3>(cuNDArray< vector_td<float,3> > *);
template  thrust::device_vector<float> Gadgetron::test_max<float,4>(cuNDArray< vector_td<float,4> > *);

template  thrust::device_vector<double> Gadgetron::test_max<double,1>(cuNDArray< vector_td<double,1> > *);
template thrust::device_vector<double> Gadgetron::test_max<double,2>(cuNDArray< vector_td<double,2> > *);
template thrust::device_vector<double> Gadgetron::test_max<double,3>(cuNDArray< vector_td<double,3> > *);
template thrust::device_vector<double> Gadgetron::test_max<double,4>(cuNDArray< vector_td<double,4> > *);
