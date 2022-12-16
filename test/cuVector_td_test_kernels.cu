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
	cudaDeviceSynchronize();
	CHECK_FOR_CUDA_ERROR();
}


template<typename T, unsigned int D>
struct test_norm_functor : public thrust::unary_function<T,vector_td<T,D> >
{
 __host__ __device__ T operator()(const vector_td<T,D> &x) const {return norm(x);}
};
template<class T, unsigned int D> thrust::device_vector<T> Gadgetron::test_norm(cuNDArray< vector_td<T,D> >* data){

	thrust::device_vector<T> out(data->get_number_of_elements());
	thrust::transform(data->begin(),data->end(),out.begin(),test_norm_functor<T,D>());
	cudaDeviceSynchronize();
	CHECK_FOR_CUDA_ERROR();
	return out;
}



template<typename T, unsigned int D>
struct test_min_functor : public thrust::unary_function<T,vector_td<T,D> >
{
 __host__ __device__ T operator()(const vector_td<T,D> &x) const {return min(x);}
};
template<class T, unsigned int D> thrust::device_vector<T> Gadgetron::test_min(cuNDArray< vector_td<T,D> >* data){

	thrust::device_vector<T> out(data->get_number_of_elements());
	thrust::transform(data->begin(),data->end(),out.begin(),test_min_functor<T,D>());
	cudaDeviceSynchronize();
	CHECK_FOR_CUDA_ERROR();
	return out;
}


template<typename T, unsigned int D>
struct test_max_functor : public thrust::unary_function<T,vector_td<T,D> >
{
 __host__ __device__ T operator()(const vector_td<T,D> &x) const {return max(x);}
};
template<class T, unsigned int D> thrust::device_vector<T> Gadgetron::test_max(cuNDArray< vector_td<T,D> >* data){

	thrust::device_vector<T> out(data->get_number_of_elements());
	thrust::transform(data->begin(),data->end(),out.begin(),test_max_functor<T,D>());
	cudaDeviceSynchronize();
	CHECK_FOR_CUDA_ERROR();
	return out;
}

template<typename T, unsigned int D>
struct test_amin_functor : public thrust::binary_function<vector_td<T,D>, vector_td<T,D>, vector_td<T,D> >
{
	__host__ __device__ vector_td<T,D> operator()(const vector_td<T,D> &x, const vector_td<T,D> &y) const {return amin(x,y);}

};

template<class T, unsigned int D> boost::shared_ptr<cuNDArray<vector_td<T,D> > > Gadgetron::test_amin(cuNDArray< vector_td<T,D> >* data1, cuNDArray< vector_td<T,D> >* data2){
	boost::shared_ptr<cuNDArray<vector_td<T,D> > > out( new cuNDArray<vector_td<T,D> >(data1->get_dimensions()));
	thrust::transform(data1->begin(),data1->end(),data2->begin(),out->begin(),test_amin_functor<T,D>());
	return out;
}


template<typename T, unsigned int D>
struct test_amax_functor : public thrust::binary_function<vector_td<T,D>, vector_td<T,D>, vector_td<T,D> >
{
	__host__ __device__ vector_td<T,D> operator()(const vector_td<T,D> &x, const vector_td<T,D> &y) const {return amax(x,y);}

};

template<class T, unsigned int D> boost::shared_ptr<cuNDArray<vector_td<T,D> > > Gadgetron::test_amax(cuNDArray< vector_td<T,D> >* data1, cuNDArray< vector_td<T,D> >* data2){
	boost::shared_ptr<cuNDArray<vector_td<T,D> > > out( new cuNDArray<vector_td<T,D> >(data1->get_dimensions()));
	thrust::transform(data1->begin(),data1->end(),data2->begin(),out->begin(),test_amax_functor<T,D>());
	return out;
}

template<typename T, unsigned int D>
class test_amin2_functor : public thrust::unary_function<vector_td<T,D>, vector_td<T,D> >
{
public:
	test_amin2_functor(T _val): val(_val){};
	__host__ __device__ vector_td<T,D> operator()(const vector_td<T,D> &x) const {return amin(x,val);}
	T val;
};

template<class T, unsigned int D> boost::shared_ptr<cuNDArray<vector_td<T,D> > > Gadgetron::test_amin2(cuNDArray< vector_td<T,D> >* data1, T val){
	boost::shared_ptr<cuNDArray<vector_td<T,D> > > out( new cuNDArray<vector_td<T,D> >(data1->get_dimensions()));
	thrust::transform(data1->begin(),data1->end(),out->begin(),test_amin2_functor<T,D>(val));
	return out;
}


template<typename T, unsigned int D>
class test_amax2_functor : public thrust::unary_function<vector_td<T,D>, vector_td<T,D> >
{
public:
	test_amax2_functor(T _val): val(_val){};
	__host__ __device__ vector_td<T,D> operator()(const vector_td<T,D> &x) const {return amax(x,val);}
	T val;
};

template<class T, unsigned int D> boost::shared_ptr<cuNDArray<vector_td<T,D> > > Gadgetron::test_amax2(cuNDArray< vector_td<T,D> >* data1, T val){
	boost::shared_ptr<cuNDArray<vector_td<T,D> > > out( new cuNDArray<vector_td<T,D> >(data1->get_dimensions()));
	thrust::transform(data1->begin(),data1->end(),out->begin(),test_amax2_functor<T,D>(val));
	return out;
}



template<class T, unsigned int D> void Gadgetron::vector_fill(cuNDArray< vector_td<T,D> >* data,  vector_td<T,D> val){
	thrust::fill(data->begin(),data->end(),val);
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



template boost::shared_ptr<cuNDArray<vector_td<float,1> > > Gadgetron::test_amin<float,1>(cuNDArray< vector_td<float,1> > *,cuNDArray< vector_td<float,1> > *);
template boost::shared_ptr<cuNDArray<vector_td<float,2> > > Gadgetron::test_amin<float,2>(cuNDArray< vector_td<float,2> > *, cuNDArray< vector_td<float,2> > *);
template  boost::shared_ptr<cuNDArray<vector_td<float,3> > > Gadgetron::test_amin<float,3>(cuNDArray< vector_td<float,3> > *, cuNDArray< vector_td<float,3> > *);
template  boost::shared_ptr<cuNDArray<vector_td<float,4> > > Gadgetron::test_amin<float,4>(cuNDArray< vector_td<float,4> > *, cuNDArray< vector_td<float,4> > *);

template  boost::shared_ptr<cuNDArray<vector_td<double,1> > > Gadgetron::test_amin<double,1>(cuNDArray< vector_td<double,1> > *, cuNDArray< vector_td<double,1> > *);
template boost::shared_ptr<cuNDArray<vector_td<double,2> > > Gadgetron::test_amin<double,2>(cuNDArray< vector_td<double,2> > *, cuNDArray< vector_td<double,2> > *);
template boost::shared_ptr<cuNDArray<vector_td<double,3> > > Gadgetron::test_amin<double,3>(cuNDArray< vector_td<double,3> > *, cuNDArray< vector_td<double,3> > *);
template boost::shared_ptr<cuNDArray<vector_td<double,4> > > Gadgetron::test_amin<double,4>(cuNDArray< vector_td<double,4> > *, cuNDArray< vector_td<double,4> > *);



template boost::shared_ptr<cuNDArray<vector_td<float,1> > > Gadgetron::test_amin2<float,1>(cuNDArray< vector_td<float,1> > *, float );
template boost::shared_ptr<cuNDArray<vector_td<float,2> > > Gadgetron::test_amin2<float,2>(cuNDArray< vector_td<float,2> > *, float);
template  boost::shared_ptr<cuNDArray<vector_td<float,3> > > Gadgetron::test_amin2<float,3>(cuNDArray< vector_td<float,3> > *, float);
template  boost::shared_ptr<cuNDArray<vector_td<float,4> > > Gadgetron::test_amin2<float,4>(cuNDArray< vector_td<float,4> > *, float);

template  boost::shared_ptr<cuNDArray<vector_td<double,1> > > Gadgetron::test_amin2<double,1>(cuNDArray< vector_td<double,1> > *, double);
template boost::shared_ptr<cuNDArray<vector_td<double,2> > > Gadgetron::test_amin2<double,2>(cuNDArray< vector_td<double,2> > *, double);
template boost::shared_ptr<cuNDArray<vector_td<double,3> > > Gadgetron::test_amin2<double,3>(cuNDArray< vector_td<double,3> > *, double);
template boost::shared_ptr<cuNDArray<vector_td<double,4> > > Gadgetron::test_amin2<double,4>(cuNDArray< vector_td<double,4> > *, double);



template boost::shared_ptr<cuNDArray<vector_td<float,1> > > Gadgetron::test_amax<float,1>(cuNDArray< vector_td<float,1> > *,cuNDArray< vector_td<float,1> > *);
template boost::shared_ptr<cuNDArray<vector_td<float,2> > > Gadgetron::test_amax<float,2>(cuNDArray< vector_td<float,2> > *, cuNDArray< vector_td<float,2> > *);
template  boost::shared_ptr<cuNDArray<vector_td<float,3> > > Gadgetron::test_amax<float,3>(cuNDArray< vector_td<float,3> > *, cuNDArray< vector_td<float,3> > *);
template  boost::shared_ptr<cuNDArray<vector_td<float,4> > > Gadgetron::test_amax<float,4>(cuNDArray< vector_td<float,4> > *, cuNDArray< vector_td<float,4> > *);

template  boost::shared_ptr<cuNDArray<vector_td<double,1> > > Gadgetron::test_amax<double,1>(cuNDArray< vector_td<double,1> > *, cuNDArray< vector_td<double,1> > *);
template boost::shared_ptr<cuNDArray<vector_td<double,2> > > Gadgetron::test_amax<double,2>(cuNDArray< vector_td<double,2> > *, cuNDArray< vector_td<double,2> > *);
template boost::shared_ptr<cuNDArray<vector_td<double,3> > > Gadgetron::test_amax<double,3>(cuNDArray< vector_td<double,3> > *, cuNDArray< vector_td<double,3> > *);
template boost::shared_ptr<cuNDArray<vector_td<double,4> > > Gadgetron::test_amax<double,4>(cuNDArray< vector_td<double,4> > *, cuNDArray< vector_td<double,4> > *);


template boost::shared_ptr<cuNDArray<vector_td<float,1> > > Gadgetron::test_amax2<float,1>(cuNDArray< vector_td<float,1> > *, float );
template boost::shared_ptr<cuNDArray<vector_td<float,2> > > Gadgetron::test_amax2<float,2>(cuNDArray< vector_td<float,2> > *, float);
template  boost::shared_ptr<cuNDArray<vector_td<float,3> > > Gadgetron::test_amax2<float,3>(cuNDArray< vector_td<float,3> > *, float);
template  boost::shared_ptr<cuNDArray<vector_td<float,4> > > Gadgetron::test_amax2<float,4>(cuNDArray< vector_td<float,4> > *, float);

template  boost::shared_ptr<cuNDArray<vector_td<double,1> > > Gadgetron::test_amax2<double,1>(cuNDArray< vector_td<double,1> > *, double);
template boost::shared_ptr<cuNDArray<vector_td<double,2> > > Gadgetron::test_amax2<double,2>(cuNDArray< vector_td<double,2> > *, double);
template boost::shared_ptr<cuNDArray<vector_td<double,3> > > Gadgetron::test_amax2<double,3>(cuNDArray< vector_td<double,3> > *, double);
template boost::shared_ptr<cuNDArray<vector_td<double,4> > > Gadgetron::test_amax2<double,4>(cuNDArray< vector_td<double,4> > *, double);



template void Gadgetron::vector_fill<float,1>(cuNDArray< vector_td<float,1> > *, vector_td<float,1>);
template void Gadgetron::vector_fill<float,2>(cuNDArray< vector_td<float,2> > *, vector_td<float,2>);
template void Gadgetron::vector_fill<float,3>(cuNDArray< vector_td<float,3> > *, vector_td<float,3>);
template void Gadgetron::vector_fill<float,4>(cuNDArray< vector_td<float,4> > *, vector_td<float,4>);


template void Gadgetron::vector_fill<double,1>(cuNDArray< vector_td<double,1> > *, vector_td<double,1>);
template void Gadgetron::vector_fill<double,2>(cuNDArray< vector_td<double,2> > *, vector_td<double,2>);
template void Gadgetron::vector_fill<double,3>(cuNDArray< vector_td<double,3> > *, vector_td<double,3>);
template void Gadgetron::vector_fill<double,4>(cuNDArray< vector_td<double,4> > *, vector_td<double,4>);
