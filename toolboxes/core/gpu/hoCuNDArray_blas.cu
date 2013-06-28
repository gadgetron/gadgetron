#include "hoCuNDArray_blas.h"
#include "cuNDArray_blas.h"
#include "complext.h"
#include "check_CUDA.h"

using namespace Gadgetron;

#define CUBLAS_CALL(fun) {cublasStatus_t err = fun; if (err != CUBLAS_STATUS_SUCCESS) {throw cuda_error(getCublasErrorString(err));}}

template<class T> EXPORTGPUCORE void Gadgetron::axpy( T a, hoCuNDArray<T>* x, hoCuNDArray<T>* y )
{
  int device = cudaDeviceManager::Instance()->getCurrentDevice();
  size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
  size_t batchSize = 1024*1024*(free/(sizeof(T)*2*1024*1024)); //Ensure 1Mb allocations
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

    CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(cuY.get_data_ptr(),y_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

    CUBLAS_CALL(cublas_axpy(cudaDeviceManager::Instance()->getHandle(device), curSize,
			    &a, cuX.get_data_ptr(), 1, cuY.get_data_ptr(), 1));

    CUDA_CALL(cudaMemcpy(y_ptr,cuY.get_data_ptr(),curSize*sizeof(T),cudaMemcpyDeviceToHost));
    remaining -= batchSize;
  }
}

template<class T> EXPORTGPUCORE void Gadgetron::axpy( T a, hoCuNDArray< complext<T> >*x, hoCuNDArray< complext<T> > *y )
{
  axpy( complext<T>(a), x, y );
}

template<class T> EXPORTGPUCORE T Gadgetron::dot( hoCuNDArray<T> *x, hoCuNDArray<T> *y, bool cc )
{
  int device = cudaDeviceManager::Instance()->getCurrentDevice();
  size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
  size_t batchSize = 1024*1024*(free/(sizeof(T)*2*1024*1024)); //Ensure 1Mb allocations
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

    CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(cuY.get_data_ptr(),y_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

    T cur_ret;
    CUBLAS_CALL(cublas_dot( cudaDeviceManager::Instance()->getHandle(device), curSize,
			    cuX.get_data_ptr(), 1,
			    cuY.get_data_ptr(), 1,
			    &cur_ret, cc ));

    remaining -= batchSize;
    ret += cur_ret;
  }
  return ret;
}

template<class T> EXPORTGPUCORE typename realType<T>::Type Gadgetron::nrm2( hoCuNDArray<T>* x )
{
  typedef typename realType<T>::Type REAL;
  int device = cudaDeviceManager::Instance()->getCurrentDevice();
  size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
  size_t batchSize = 1024*1024*(free/(sizeof(T)*1024*1024)); //Ensure 1Mb allocations
  size_t remaining = x->get_number_of_elements();
  batchSize = std::min(batchSize,remaining);
  T* x_ptr = x->get_data_ptr();
  std::vector<unsigned int> dims;
  dims.push_back(batchSize);
  cuNDArray<T> cuX(&dims);
  REAL ret = 0;

  for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

    size_t curSize = std::min(batchSize,remaining);
    CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

    REAL cur_ret;
    CUBLAS_CALL(cublas_nrm2<T>( cudaDeviceManager::Instance()->getHandle(device), batchSize,
				cuX.get_data_ptr(), 1, &cur_ret));
    remaining -= batchSize;
    ret += cur_ret*cur_ret;
  }
  return std::sqrt(ret);
}

template<class T> EXPORTGPUCORE typename realType<T>::Type Gadgetron::asum( hoCuNDArray<T>* x )
{
  typedef typename realType<T>::Type REAL;
  int device = cudaDeviceManager::Instance()->getCurrentDevice();
  size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
  size_t batchSize = 1024*1024*(free/(sizeof(T)*1024*1024)); //Ensure 1Mb allocations
  size_t remaining = x->get_number_of_elements();
  batchSize = std::min(batchSize,remaining);
  T* x_ptr = x->get_data_ptr();
  std::vector<unsigned int> dims;
  dims.push_back(batchSize);
  cuNDArray<T> cuX(&dims);
  REAL ret = 0;

  for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

    size_t curSize = std::min(batchSize,remaining);
    CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

    REAL cur_ret;
    CUBLAS_CALL(cublas_asum( cudaDeviceManager::Instance()->getHandle(device), batchSize,
			     cuX.get_data_ptr(), 1,
			     &cur_ret));
    remaining -= batchSize;
    ret += cur_ret;
  }
  return ret;
}

template<class T> EXPORTGPUCORE int Gadgetron::amin( hoCuNDArray<T>* x )
{
  int device = cudaDeviceManager::Instance()->getCurrentDevice();
  size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
  size_t batchSize = 1024*1024*(free/(sizeof(T)*1024*1024)); //Ensure 1Mb allocations
  size_t remaining = x->get_number_of_elements();
  batchSize = std::min(batchSize,remaining);
  T* x_ptr = x->get_data_ptr();
  std::vector<unsigned int> dims;
  dims.push_back(batchSize);
  cuNDArray<T> cuX(&dims);
  std::vector<int> results;
 
 for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

    size_t curSize = std::min(batchSize,remaining);
    CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

    int cur_ret;
    CUBLAS_CALL(cublas_amin( cudaDeviceManager::Instance()->getHandle(device), batchSize,
			     cuX.get_data_ptr(), 1,
			     &cur_ret));
    remaining -= batchSize;
    results.push_back(cur_ret+i*batchSize-1);
  }

  int res =0;
  for (int i =0; i < results.size(); i++){
    if (abs(x_ptr[results[i]]) < abs(x_ptr[res])) res = results[i];
  }
  return res;
}

template<class T> EXPORTGPUCORE int Gadgetron::amax( hoCuNDArray<T>* x )
{
  int device = cudaDeviceManager::Instance()->getCurrentDevice();
  size_t free = cudaDeviceManager::Instance()->getFreeMemory(device);
  size_t batchSize = 1024*1024*(free/(sizeof(T)*1024*1024)); //Ensure 1Mb allocations
  size_t remaining = x->get_number_of_elements();
  batchSize = std::min(batchSize,remaining);
  T* x_ptr = x->get_data_ptr();
  std::vector<unsigned int> dims;
  dims.push_back(batchSize);
  cuNDArray<T> cuX(&dims);
  std::vector<int> results;

  for (size_t i = 0; i < (x->get_number_of_elements()-1)/batchSize+1; i++){

    size_t curSize = std::min(batchSize,remaining);
    CUDA_CALL(cudaMemcpy(cuX.get_data_ptr(),x_ptr+i*batchSize,curSize*sizeof(T),cudaMemcpyHostToDevice));

    int cur_ret;
    CUBLAS_CALL(cublas_amax( cudaDeviceManager::Instance()->getHandle(device), batchSize,
			     cuX.get_data_ptr(), 1,
			     &cur_ret));
    remaining -= batchSize;
    results.push_back(cur_ret+i*batchSize-1);
  }

  int res =0;
  for (int i =0; i < results.size(); i++){
    if (abs(x_ptr[results[i]]) > abs(x_ptr[res])) res = results[i];
  }
  return res;
}

//
// Instantiation
//

template EXPORTGPUCORE float Gadgetron::dot(hoCuNDArray<float>*,hoCuNDArray<float>*,bool);
template EXPORTGPUCORE float Gadgetron::nrm2(hoCuNDArray<float>*);
template EXPORTGPUCORE void Gadgetron::axpy(float,hoCuNDArray<float>*,hoCuNDArray<float>*);
template EXPORTGPUCORE int Gadgetron::amin(hoCuNDArray<float>*);
template EXPORTGPUCORE int Gadgetron::amax(hoCuNDArray<float>*);
template EXPORTGPUCORE float Gadgetron::asum(hoCuNDArray<float>*);

template EXPORTGPUCORE double Gadgetron::dot(hoCuNDArray<double>*,hoCuNDArray<double>*,bool);
template EXPORTGPUCORE double Gadgetron::nrm2(hoCuNDArray<double>*);
template EXPORTGPUCORE void Gadgetron::axpy(double,hoCuNDArray<double>*,hoCuNDArray<double>*);
template EXPORTGPUCORE int Gadgetron::amin(hoCuNDArray<double>*);
template EXPORTGPUCORE int Gadgetron::amax(hoCuNDArray<double>*);
template EXPORTGPUCORE double Gadgetron::asum(hoCuNDArray<double>*);

template EXPORTGPUCORE float_complext Gadgetron::dot(hoCuNDArray<float_complext>*,hoCuNDArray<float_complext>*,bool);
template EXPORTGPUCORE float Gadgetron::nrm2(hoCuNDArray<float_complext>*);
template EXPORTGPUCORE void Gadgetron::axpy(float_complext,hoCuNDArray<float_complext>*,hoCuNDArray<float_complext>*);
template EXPORTGPUCORE void Gadgetron::axpy(float,hoCuNDArray<float_complext>*,hoCuNDArray<float_complext>*);
template EXPORTGPUCORE int Gadgetron::amin(hoCuNDArray<float_complext>*);
template EXPORTGPUCORE int Gadgetron::amax(hoCuNDArray<float_complext>*);
template EXPORTGPUCORE float Gadgetron::asum(hoCuNDArray<float_complext>*);

template EXPORTGPUCORE double_complext Gadgetron::dot(hoCuNDArray<double_complext>*,hoCuNDArray<double_complext>*,bool);
template EXPORTGPUCORE double Gadgetron::nrm2(hoCuNDArray<double_complext>*);
template EXPORTGPUCORE void Gadgetron::axpy(double_complext,hoCuNDArray<double_complext>*,hoCuNDArray<double_complext>*);
template EXPORTGPUCORE void Gadgetron::axpy(double,hoCuNDArray<double_complext>*,hoCuNDArray<double_complext>*);
template EXPORTGPUCORE int Gadgetron::amin(hoCuNDArray<double_complext>*);
template EXPORTGPUCORE int Gadgetron::amax(hoCuNDArray<double_complext>*);
template EXPORTGPUCORE double Gadgetron::asum(hoCuNDArray<double_complext>*);
