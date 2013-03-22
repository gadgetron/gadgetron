#include "hoCuGTBLAS.h"

#include "cuGTBLAS.h"

#include "complext.h"
#include "check_CUDA.h"

using namespace Gadgetron;

#define CUBLAS_CALL(fun) {cublasStatus_t err = fun; if (err != CUBLAS_STATUS_SUCCESS) {BOOST_THROW_EXCEPTION(cuda_error(getCublasErrorString(err)));}}

template<class T> void Gadgetron::axpy(T a, hoCuNDArray<T>* x, hoCuNDArray<T>* y, int device){


		size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
		size_t batchSize = (free/(sizeof(T)*2));

		size_t remaining = x->get_number_of_elements();
		batchSize = std::min(batchSize,remaining);
		T* x_ptr = x->get_data_ptr();
		T* y_ptr = y->get_data_ptr();
		std::vector<unsigned int> dims;
		dims.push_back(batchSize);
		cuNDArray<T> cuX(&dims);
		cuNDArray<T> cuY(&dims);

		for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

			size_t curSize = std::min(batchSize,remaining);


			CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));
			CUDA_CALL(cudaMemcpy(cuY.get_data_ptr(),y_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

		  CUBLAS_CALL(cublas_axpy(cudaDeviceManager::Instance()->getHandle(device), curSize,
				  &a, cuX.get_data_ptr(), 1,
				  cuY.get_data_ptr(), 1));
			CUDA_CALL(cudaMemcpy(y_ptr,cuY.get_data_ptr(),curSize*sizeof(T),cudaMemcpyDeviceToHost));

			remaining -= batchSize;

		}
}

template<class T> T Gadgetron::dot( hoCuNDArray<T>* x, hoCuNDArray<T>* y, int device){
		size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
		size_t batchSize = (free/(sizeof(T)*2));

		size_t remaining = x->get_number_of_elements();
		batchSize = std::min(batchSize,remaining);
		T* x_ptr = x->get_data_ptr();
		T* y_ptr = y->get_data_ptr();
		std::vector<unsigned int> dims;
		dims.push_back(batchSize);
		cuNDArray<T> cuX(&dims);
		cuNDArray<T> cuY(&dims);
		T ret = T(0);
		for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

			size_t curSize = std::min(batchSize,remaining);


			CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));
			CUDA_CALL(cudaMemcpy(cuY.get_data_ptr(),y_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

			T cur_ret;

			  CUBLAS_CALL(cublas_dot( cudaDeviceManager::Instance()->getHandle(device), curSize,
					   cuX.get_data_ptr(), 1,
					   cuY.get_data_ptr(), 1,
					   &cur_ret));

			remaining -= batchSize;
			ret += cur_ret;
		}
		return ret;
}



template<class T> typename realType<T>::type Gadgetron::nrm2( hoCuNDArray<T>* x, int device){
		typedef typename realType<T>::type REAL;
		size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
		size_t batchSize = (free/(sizeof(T)));

		size_t remaining = x->get_number_of_elements();
		batchSize = std::min(batchSize,remaining);
		T* x_ptr = x->get_data_ptr();
		std::vector<unsigned int> dims;
		dims.push_back(batchSize);
		cuNDArray<T> cuX(&dims);
		REAL ret = 0;
		for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

			size_t curSize = std::min(batchSize,remaining);
			CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

			REAL cur_ret;

			CUBLAS_CALL(cublas_nrm2<T>( cudaDeviceManager::Instance()->getHandle(device), batchSize,
					    cuX.get_data_ptr(), 1,
					   &cur_ret));
			remaining -= batchSize;
			ret += cur_ret*cur_ret;
		}
		return std::sqrt(ret);
}


template<class T> typename realType<T>::type Gadgetron::asum(hoCuNDArray<T>* x,int device){
		typedef typename realType<T>::type REAL;
		size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
		size_t batchSize = (free/(sizeof(T)));

		size_t remaining = x->get_number_of_elements();
		batchSize = std::min(batchSize,remaining);
		T* x_ptr = x->get_data_ptr();
		std::vector<unsigned int> dims;
		dims.push_back(batchSize);
		cuNDArray<T> cuX(&dims);
		REAL ret = 0;
		for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

			size_t curSize = std::min(batchSize,remaining);
			CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

			REAL cur_ret;

			CUBLAS_CALL(cublas_asum( cudaDeviceManager::Instance()->getHandle(device), batchSize,
					    cuX.get_data_ptr(), 1,
					   &cur_ret));
			remaining -= batchSize;
			ret += cur_ret;
		}
		return ret;
}

template<class T> int Gadgetron::amin(hoCuNDArray<T>* x,int device){

		size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
		size_t batchSize = (free/(sizeof(T)));

		size_t remaining = x->get_number_of_elements();
		batchSize = std::min(batchSize,remaining);
		T* x_ptr = x->get_data_ptr();
		std::vector<unsigned int> dims;
		dims.push_back(batchSize);
		cuNDArray<T> cuX(&dims);
		std::vector<int> results;
		for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

			size_t curSize = std::min(batchSize,remaining);
			CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

			int cur_ret;

			CUBLAS_CALL(cublas_amin( cudaDeviceManager::Instance()->getHandle(device), batchSize,
					    cuX.get_data_ptr(), 1,
					   &cur_ret));
			remaining -= batchSize;
			results.push_back(cur_ret+i*curSize-1);
		}
		int res =0;
		for (int i =0; i < results.size(); i++){
			if (abs(x_ptr[results[i]]) < abs(x_ptr[res])) res = results[i];
		}
		return res;

}

template<class T> int Gadgetron::amax(hoCuNDArray<T>* x,int device){

		size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
		size_t batchSize = (free/(sizeof(T)));

		size_t remaining = x->get_number_of_elements();
		batchSize = std::min(batchSize,remaining);
		T* x_ptr = x->get_data_ptr();
		std::vector<unsigned int> dims;
		dims.push_back(batchSize);
		cuNDArray<T> cuX(&dims);
		std::vector<int> results;
		for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

			size_t curSize = std::min(batchSize,remaining);
			CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*curSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

			int cur_ret;

			CUBLAS_CALL(cublas_amax( cudaDeviceManager::Instance()->getHandle(device), batchSize,
					    cuX.get_data_ptr(), 1,
					   &cur_ret));
			remaining -= batchSize;
			results.push_back(cur_ret+i*curSize-1);
		}
		int res =0;
		for (int i =0; i < results.size(); i++){
			if (abs(x_ptr[results[i]]) > abs(x_ptr[res])) res = results[i];
		}
		return res;

}
template float Gadgetron::dot(hoCuNDArray<float> *x,hoCuNDArray<float> *y,int device);
template float Gadgetron::nrm2( hoCuNDArray<float>* arr, int device);
template void Gadgetron::axpy(float a, hoCuNDArray<float>* x, hoCuNDArray<float>* y,int device);
template int Gadgetron::amin(hoCuNDArray<float>* x,int device);
template int Gadgetron::amax(hoCuNDArray<float>* x,int device);
template float Gadgetron::asum(hoCuNDArray<float>* x,int device);

template double Gadgetron::dot(hoCuNDArray<double> *x,hoCuNDArray<double> *y,int device);
template double Gadgetron::nrm2( hoCuNDArray<double>* arr, int device);
template void Gadgetron::axpy(double a, hoCuNDArray<double>* x, hoCuNDArray<double>* y,int device);
template int Gadgetron::amin(hoCuNDArray<double>* x,int device);
template int Gadgetron::amax(hoCuNDArray<double>* x,int device);
template double Gadgetron::asum(hoCuNDArray<double>* x,int device);

template float_complext Gadgetron::dot(hoCuNDArray<float_complext> *x,hoCuNDArray<float_complext> *y,int device);
template float Gadgetron::nrm2( hoCuNDArray<float_complext>* arr, int device);
template void Gadgetron::axpy(float_complext a, hoCuNDArray<float_complext>* x, hoCuNDArray<float_complext>* y,int device);
template void Gadgetron::axpy(float a, hoCuNDArray<float_complext>* x, hoCuNDArray<float_complext>* y,int device);

template int Gadgetron::amin(hoCuNDArray<float_complext>* x,int device);
template int Gadgetron::amax(hoCuNDArray<float_complext>* x,int device);
template float Gadgetron::asum(hoCuNDArray<float_complext>* x,int device);


template double_complext Gadgetron::dot(hoCuNDArray<double_complext> *x,hoCuNDArray<double_complext> *y,int device);
template double Gadgetron::nrm2( hoCuNDArray<double_complext>* arr, int device);
template void Gadgetron::axpy(double_complext a, hoCuNDArray<double_complext>* x, hoCuNDArray<double_complext>* y,int device);
template void Gadgetron::axpy(double a, hoCuNDArray<double_complext>* x, hoCuNDArray<double_complext>* y,int device);

template int Gadgetron::amin(hoCuNDArray<double_complext>* x,int device);
template int Gadgetron::amax(hoCuNDArray<double_complext>* x,int device);
template double Gadgetron::asum(hoCuNDArray<double_complext>* x,int device);
