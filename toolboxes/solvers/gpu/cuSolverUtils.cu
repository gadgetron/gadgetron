#include "complext.h"
#include "cuSolverUtils.h"
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include "cuNDArray_math.h"
#define MAX_THREADS_PER_BLOCK 512

using namespace Gadgetron;
template <class T> __global__ static void filter_kernel(T* x, T* g, int elements){
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < elements){
		if ( x[idx] <= T(0) && g[idx] > 0) g[idx]=T(0);
	}
}

template <class REAL> __global__ static void filter_kernel(complext<REAL>* x, complext<REAL>* g, int elements){
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < elements){
		if ( real(x[idx]) <= REAL(0) && real(g[idx]) > 0) g[idx]._real = REAL(0);
		g[idx]._imag=REAL(0);
	}
}

template <class T> void EXPORTGPUSOLVERS Gadgetron::solver_non_negativity_filter(cuNDArray<T>* x , cuNDArray<T>* g)
{
	int elements = g->get_number_of_elements();

	int threadsPerBlock = std::min(elements,MAX_THREADS_PER_BLOCK);
	dim3 dimBlock( threadsPerBlock);
	int totalBlocksPerGrid = std::max(1,elements/MAX_THREADS_PER_BLOCK);
	dim3 dimGrid(totalBlocksPerGrid);

	filter_kernel<typename realType<T>::Type><<<dimGrid,dimBlock>>>(x->get_data_ptr(),g->get_data_ptr(),elements);
}



template<class T> struct updateF_functor{

	typedef typename realType<T>::Type REAL;
	updateF_functor(REAL alpha_, REAL sigma_){

		alpha= alpha_;
		sigma = sigma_;
	}
	__device__ __inline__ T operator() (T val){
		return val/(1+alpha*sigma)/max(REAL(1),abs(val/(1+alpha*sigma)));
	}
	typename realType<T>::Type alpha, sigma;
};

template<class T>
inline void Gadgetron::updateF(cuNDArray<T>& data,
		typename realType<T>::Type alpha, typename realType<T>::Type sigma) {
	thrust::transform(data.begin(),data.end(),data.begin(),updateF_functor<T>(alpha,sigma));
}


template<class T> struct updateFgroup_functor {

	typedef typename realType<T>::Type REAL;
	updateFgroup_functor(REAL alpha_, REAL sigma_) : alpha(alpha_), sigma(sigma_){
	}

	__device__ __inline__ T operator() (thrust::tuple<T,typename realType<T>::Type> tup){
			return thrust::get<0>(tup)/(1+alpha*sigma)/max(REAL(1),thrust::get<1>(tup)/(1+alpha*sigma));
	}

	typename realType<T>::Type alpha, sigma;
};

template<class T> struct add_square_functor{

	__device__ __inline__ typename realType<T>::Type operator() (thrust::tuple<T,typename realType<T>::Type> tup){
		T val = thrust::get<0>(tup);
		return thrust::get<1>(tup)+norm(val);
	}
};
template<class T>
inline void Gadgetron::updateFgroup(std::vector<cuNDArray<T> >& datas,
		typename realType<T>::Type alpha, typename realType<T>::Type sigma) {

	cuNDArray<typename realType<T>::Type> squares(datas.front().get_dimensions());
	clear(&squares);
	for (int i = 0; i < datas.size(); i++)
		thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(datas[i].begin(),squares.begin())),
				thrust::make_zip_iterator(thrust::make_tuple(datas[i].end(),squares.end())), squares.begin(), add_square_functor<T>());

	sqrt_inplace(&squares);
	for (int i = 0 ; i < datas.size(); i++){
		thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(datas[i].begin(),squares.begin())),
				thrust::make_zip_iterator(thrust::make_tuple(datas[i].end(),squares.end())), datas[i].begin(), updateFgroup_functor<T>(alpha,sigma));
	}
}

template void EXPORTGPUSOLVERS Gadgetron::updateF<float>(cuNDArray<float>& data, float alpha, float sigma);
template void EXPORTGPUSOLVERS Gadgetron::updateF<double>(cuNDArray<double>& data,double alpha, double sigma);
template void EXPORTGPUSOLVERS Gadgetron::updateF<float_complext>(cuNDArray<float_complext>& data, float alpha, float sigma);
template void EXPORTGPUSOLVERS Gadgetron::updateF<double_complext>(cuNDArray<double_complext>& data, double alpha, double sigma);

template void EXPORTGPUSOLVERS Gadgetron::updateFgroup<float>(std::vector<cuNDArray<float> >& data, float alpha, float sigma);
template void EXPORTGPUSOLVERS Gadgetron::updateFgroup<double>(std::vector<cuNDArray<double> >& data,double alpha, double sigma);
template void EXPORTGPUSOLVERS Gadgetron::updateFgroup<float_complext>(std::vector<cuNDArray<float_complext> >& data, float alpha, float sigma);
template void EXPORTGPUSOLVERS Gadgetron::updateFgroup<double_complext>(std::vector<cuNDArray<double_complext> >& data, double alpha, double sigma);


template void EXPORTGPUSOLVERS Gadgetron::solver_non_negativity_filter<float>(cuNDArray<float>*, cuNDArray<float>*);
template void EXPORTGPUSOLVERS Gadgetron::solver_non_negativity_filter<double>(cuNDArray<double>*, cuNDArray<double>*);
template void EXPORTGPUSOLVERS Gadgetron::solver_non_negativity_filter<float_complext>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);
template void EXPORTGPUSOLVERS Gadgetron::solver_non_negativity_filter<double_complext>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);

