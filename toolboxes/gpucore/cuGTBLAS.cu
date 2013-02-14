#include <cublas_v2.h>
#include "cuGTBLAS.h"
#include "GadgetronCuException.h"
#include "complext.h"
#include "cudaDeviceManager.h"

namespace Gadgetron{

#define CUBLAS_CALL(fun) {cublasStatus_t err = fun; if (err != CUBLAS_STATUS_SUCCESS) {BOOST_THROW_EXCEPTION(cuda_error(getCublasErrorString(err)));}}
//NRM2

template<> cublasStatus_t cublas_nrm2<float>(cublasHandle_t hndl, int n, const float*  x, int inc, float* res){
		return cublasSnrm2(hndl,n,x,inc,res);
}
template<> cublasStatus_t cublas_nrm2<double>(cublasHandle_t hndl, int n, const double*  x, int inc, double* res){
		return cublasDnrm2(hndl,n,x,inc,res);
}
template<> cublasStatus_t cublas_nrm2<float_complext>(cublasHandle_t hndl, int n, const float_complext*  x, int inc, float* res){
		return cublasScnrm2(hndl,n,(const cuComplex*)x,inc,res);
}
template<> cublasStatus_t cublas_nrm2<double_complext>(cublasHandle_t hndl, int n, const double_complext*  x, int inc, double* res){
		return cublasDznrm2(hndl,n,(const cuDoubleComplex*) x,inc,res);
}

//DOT
template<> cublasStatus_t cublas_dot<float>(cublasHandle_t hndl, int n , const float* x , int incx, const  float* y , int incy, float* res){
	return cublasSdot( hndl, n, x, incx, y, incy, res);
}
template<> cublasStatus_t cublas_dot<double>(cublasHandle_t hndl, int n , const double* x , int incx, const  double* y , int incy, double* res){
	return cublasDdot( hndl, n, x, incx, y, incy, res);
}
template<> cublasStatus_t cublas_dot<float_complext>(cublasHandle_t hndl, int n , const float_complext* x ,
		int incx, const  float_complext* y , int incy, float_complext* res){
	return cublasCdotc( hndl, n, (const cuComplex*) x, incx, (const cuComplex*) y, incy, (cuComplex*) res);
}
template<> cublasStatus_t cublas_dot<double_complext>(cublasHandle_t hndl, int n , const double_complext* x ,
		int incx, const  double_complext* y , int incy, double_complext* res){
	return cublasZdotc( hndl, n, (const cuDoubleComplex*) x, incx, (const cuDoubleComplex*) y, incy, (cuDoubleComplex*) res);
}

// AXPY
template<> cublasStatus_t cublas_axpy<float>(cublasHandle_t hndl , int n , const float* a , const float* x , int incx ,  float* y , int incy){
	return cublasSaxpy(hndl,n,a,x,incx,y,incy);
}
template<> cublasStatus_t cublas_axpy<double>(cublasHandle_t hndl , int n , const double* a , const double* x , int incx ,  double* y , int incy){
	return cublasDaxpy(hndl,n,a,x,incx,y,incy);
}
template<> cublasStatus_t cublas_axpy<float_complext>(cublasHandle_t hndl , int n , const float_complext* a , const float_complext* x , int incx ,  float_complext* y , int incy){
	return cublasCaxpy(hndl,n,(const cuComplex*) a, (const cuComplex*) x,incx, (cuComplex*)y,incy);
}
template<> cublasStatus_t cublas_axpy<double_complext>(cublasHandle_t hndl , int n , const double_complext* a , const double_complext* x , int incx ,  double_complext* y , int incy){
	return cublasZaxpy(hndl,n,(const cuDoubleComplex*) a, (const cuDoubleComplex*) x,incx, (cuDoubleComplex*)y,incy);
}

//SUM
template<> cublasStatus_t cublas_asum<float>(cublasHandle_t hndl, int n,const float *x, int incx, float *result){
	return cublasSasum(hndl,n,x,incx,result);
}
template<> cublasStatus_t cublas_asum<double>(cublasHandle_t hndl, int n,const double *x, int incx, double *result){
	return cublasDasum(hndl,n,x,incx,result);
}
template<> cublasStatus_t cublas_asum<float_complext>(cublasHandle_t hndl, int n,const float_complext *x, int incx, float *result){
	return cublasScasum(hndl,n,(const cuComplex*) x,incx,result);
}
template<> cublasStatus_t cublas_asum<double_complext>(cublasHandle_t hndl, int n,const double_complext *x, int incx, double *result){
	return cublasDzasum(hndl,n,(const cuDoubleComplex*) x,incx,result);
}

//AMIN
template<> cublasStatus_t cublas_amin<float>(cublasHandle_t hndl, int n,const float *x, int incx, int *result){
	return cublasIsamin(hndl,n,x,incx,result);
}
template<> cublasStatus_t cublas_amin<double>(cublasHandle_t hndl, int n,const double *x, int incx, int *result){
	return cublasIdamin(hndl,n,x,incx,result);
}
template<> cublasStatus_t cublas_amin<float_complext>(cublasHandle_t hndl, int n,const float_complext *x, int incx, int *result){
	return cublasIcamin(hndl,n, (const cuComplex* ) x,incx,result);
}
template<> cublasStatus_t cublas_amin<double_complext>(cublasHandle_t hndl, int n,const double_complext *x, int incx, int *result){
	return cublasIzamin(hndl,n, (const cuDoubleComplex* ) x,incx,result);
}

//AMAX
template<> cublasStatus_t cublas_amax<float>(cublasHandle_t hndl, int n,const float *x, int incx, int *result){
	return cublasIsamax(hndl,n,x,incx,result);
}
template<> cublasStatus_t cublas_amax<double>(cublasHandle_t hndl, int n,const double *x, int incx, int *result){
	return cublasIdamax(hndl,n,x,incx,result);
}
template<> cublasStatus_t cublas_amax<float_complext>(cublasHandle_t hndl, int n,const float_complext *x, int incx, int *result){
	return cublasIcamax(hndl,n, (const cuComplex* ) x,incx,result);
}
template<> cublasStatus_t cublas_amax<double_complext>(cublasHandle_t hndl, int n,const double_complext *x, int incx, int *result){
	return cublasIzamax(hndl,n, (const cuDoubleComplex* ) x,incx,result);
}


template<class T> typename realType<T>::type
nrm2( cuNDArray<T>* arr, int device )
{
	typedef typename realType<T>::type REAL;
	REAL ret;

  CUBLAS_CALL(cublas_nrm2<T>( cudaDeviceManager::Instance()->getHandle(device), arr->get_number_of_elements(),
		    arr->get_data_ptr(), 1,
		   &ret));
  cudaThreadSynchronize();
  return ret;
}

template<class T> T
dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, int device )
{

  T ret;

  CUBLAS_CALL(cublas_dot( cudaDeviceManager::Instance()->getHandle(device), arr1->get_number_of_elements(),
		   arr1->get_data_ptr(), 1,
		   arr2->get_data_ptr(), 1,
		   &ret));

  cudaThreadSynchronize();
  return ret;
}


template<class T> void
axpy(T a,  cuNDArray<T>* x, cuNDArray<T>* y, int device )
{
  CUBLAS_CALL(cublas_axpy(cudaDeviceManager::Instance()->getHandle(device), x->get_number_of_elements(),
		  &a, x->get_data_ptr(), 1,
		  y->get_data_ptr(), 1));
  cudaThreadSynchronize();
}


template<class T> typename realType<T>::type asum(cuNDArray<T>* x,int device){

	typename realType<T>::type result;
	CUBLAS_CALL(cublas_asum(cudaDeviceManager::Instance()->getHandle(device),x->get_number_of_elements(),x->get_data_ptr(),1,&result));

	return result;
}

template<class T> int amin(cuNDArray<T>* x,int device){

	int result;
	CUBLAS_CALL(cublas_amin(cudaDeviceManager::Instance()->getHandle(device),x->get_number_of_elements(),x->get_data_ptr(),1,&result));
	return result;
}

template<class T> int amax(cuNDArray<T>* x,int device){
	int result;
	CUBLAS_CALL(cublas_amax(cudaDeviceManager::Instance()->getHandle(device),x->get_number_of_elements(),x->get_data_ptr(),1,&result));
	return result;
}


std::string getCublasErrorString(cublasStatus_t err){
	switch (err){
	case CUBLAS_STATUS_NOT_INITIALIZED:
		return "NOT INITIALIZED";
	case CUBLAS_STATUS_ALLOC_FAILED:
		return "ALLOC FAILED";
	case CUBLAS_STATUS_INVALID_VALUE:
		return "INVALID VALUE";
	case CUBLAS_STATUS_ARCH_MISMATCH:
		return "ARCH MISMATCH";
	case CUBLAS_STATUS_MAPPING_ERROR:
		return "MAPPING ERROR";
	case CUBLAS_STATUS_EXECUTION_FAILED:
		return "EXECUTION FAILED";
	case CUBLAS_STATUS_INTERNAL_ERROR:
		return "INTERNAL ERROR";

	case CUBLAS_STATUS_SUCCESS:
		return "SUCCES";
	default:
		return "UNKNOWN CUBLAS ERROR";
	}
}

template float dot(cuNDArray<float> *x,cuNDArray<float> *y,int device);
template float nrm2( cuNDArray<float>* arr, int device);
template void axpy(float a, cuNDArray<float>* x, cuNDArray<float>* y,int device);
template int amin(cuNDArray<float>* x,int device);
template int amax(cuNDArray<float>* x,int device);
template float asum(cuNDArray<float>* x,int device);

template double dot(cuNDArray<double> *x,cuNDArray<double> *y,int device);
template double nrm2( cuNDArray<double>* arr, int device);
template void axpy(double a, cuNDArray<double>* x, cuNDArray<double>* y,int device);
template int amin(cuNDArray<double>* x,int device);
template int amax(cuNDArray<double>* x,int device);
template double asum(cuNDArray<double>* x,int device);

template float_complext dot(cuNDArray<float_complext> *x,cuNDArray<float_complext> *y,int device);
template float nrm2( cuNDArray<float_complext>* arr, int device);
template void axpy(float_complext a, cuNDArray<float_complext>* x, cuNDArray<float_complext>* y,int device);
template void axpy(float a, cuNDArray<float_complext>* x, cuNDArray<float_complext>* y,int device);

template int amin(cuNDArray<float_complext>* x,int device);
template int amax(cuNDArray<float_complext>* x,int device);
template float asum(cuNDArray<float_complext>* x,int device);


template double_complext dot(cuNDArray<double_complext> *x,cuNDArray<double_complext> *y,int device);
template double nrm2( cuNDArray<double_complext>* arr, int device);
template void axpy(double_complext a, cuNDArray<double_complext>* x, cuNDArray<double_complext>* y,int device);
template void axpy(double a, cuNDArray<double_complext>* x, cuNDArray<double_complext>* y,int device);

template int amin(cuNDArray<double_complext>* x,int device);
template int amax(cuNDArray<double_complext>* x,int device);
template double asum(cuNDArray<double_complext>* x,int device);

}

