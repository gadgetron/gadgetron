#include <cublas_v2.h>
#include "cuGTBLAS.h"
#include "gadgetronException.h"
#include "complext.h"

// Some device properties we query once to eliminate runtime overhead
//
static int num_devices = 0;
static int *warp_size = 0x0;
static int *max_blockdim = 0x0;
static int *max_griddim = 0x0;
static cublasHandle_t *handle = 0x0;
// Initialize static variables
//
static void initialize_static_variables()
{
  // This function is executed only once
  if( num_devices ) return;

  if( cudaGetDeviceCount( &num_devices ) != cudaSuccess) {
    num_devices = 0;
    throw cuda_error("Error: no Cuda devices present.");
  }

  int old_device;
  if( cudaGetDevice(&old_device) != cudaSuccess ) {
    throw cuda_error("Error: unable to get device no.");
  }

  warp_size = new int[num_devices];
  max_blockdim = new int[num_devices];
  max_griddim = new int[num_devices];
  handle = new cublasHandle_t[num_devices];

  for( int device=0; device<num_devices; device++ ){

    if( cudaSetDevice(device) != cudaSuccess ) {
      throw cuda_error("Error: unable to set device no.");
    }

    cudaDeviceProp deviceProp;

    if( cudaGetDeviceProperties( &deviceProp, device ) != cudaSuccess) {
      throw cuda_error("Error: unable to determine device properties.");

    }

    warp_size[device] = deviceProp.warpSize;
    max_blockdim[device] = deviceProp.maxThreadsDim[0];
    max_griddim[device] = deviceProp.maxGridSize[0];

    if (cublasCreate(&handle[device]) != CUBLAS_STATUS_SUCCESS) {
      throw cuda_error("Error: unable to create cublas handle for device ");
    }

    cublasSetPointerMode( handle[device], CUBLAS_POINTER_MODE_HOST );

  }

  if( cudaSetDevice(old_device) != cudaSuccess ) {
    throw cuda_error("Error: unable to restore device no");
  }

}



//NRM2

template<> int cublas_nrm2<float>(cublasHandle_t hndl, int n, const float*  x, int inc, float* res){
		return cublasSnrm2(hndl,n,x,inc,res);
}
template<> int cublas_nrm2<double>(cublasHandle_t hndl, int n, const double*  x, int inc, double* res){
		return cublasDnrm2(hndl,n,x,inc,res);
}
template<> int cublas_nrm2<float_complext>(cublasHandle_t hndl, int n, const float_complext*  x, int inc, float* res){
		return cublasScnrm2(hndl,n,(const cuComplex*)x,inc,res);
}
template<> int cublas_nrm2<double_complext>(cublasHandle_t hndl, int n, const double_complext*  x, int inc, double* res){
		return cublasDznrm2(hndl,n,(const cuDoubleComplex*) x,inc,res);
}

//DOT
template<> int cublas_dot<float>(cublasHandle_t hndl, int n , const float* x , int incx, const  float* y , int incy, float* res){
	return cublasSdot( hndl, n, x, incx, y, incy, res);
}
template<> int cublas_dot<double>(cublasHandle_t hndl, int n , const double* x , int incx, const  double* y , int incy, double* res){
	return cublasDdot( hndl, n, x, incx, y, incy, res);
}
template<> int cublas_dot<float_complext>(cublasHandle_t hndl, int n , const float_complext* x ,
		int incx, const  float_complext* y , int incy, float_complext* res){
	return cublasCdotc( hndl, n, (const cuComplex*) x, incx, (const cuComplex*) y, incy, (cuComplex*) res);
}
template<> int cublas_dot<double_complext>(cublasHandle_t hndl, int n , const double_complext* x ,
		int incx, const  double_complext* y , int incy, double_complext* res){
	return cublasZdotc( hndl, n, (const cuDoubleComplex*) x, incx, (const cuDoubleComplex*) y, incy, (cuDoubleComplex*) res);
}

// AXPY
template<> int cublas_axpy<float>(cublasHandle_t hndl , int n , const float* a , const float* x , int incx ,  float* y , int incy){
	return cublasSaxpy(hndl,n,a,x,incx,y,incy);
}
template<> int cublas_axpy<double>(cublasHandle_t hndl , int n , const double* a , const double* x , int incx ,  double* y , int incy){
	return cublasDaxpy(hndl,n,a,x,incx,y,incy);
}
template<> int cublas_axpy<float_complext>(cublasHandle_t hndl , int n , const float_complext* a , const float_complext* x , int incx ,  float_complext* y , int incy){
	return cublasCaxpy(hndl,n,(const cuComplex*) a, (const cuComplex*) x,incx, (cuComplex*)y,incy);
}
template<> int cublas_axpy<double_complext>(cublasHandle_t hndl , int n , const double_complext* a , const double_complext* x , int incx ,  double_complext* y , int incy){
	return cublasZaxpy(hndl,n,(const cuDoubleComplex*) a, (const cuDoubleComplex*) x,incx, (cuDoubleComplex*)y,incy);
}

//SUM
template<> int cublas_asum<float>(cublasHandle_t hndl, int n,const float *x, int incx, float *result){
	return cublasSasum(hndl,n,x,incx,result);
}
template<> int cublas_asum<double>(cublasHandle_t hndl, int n,const double *x, int incx, double *result){
	return cublasDasum(hndl,n,x,incx,result);
}
template<> int cublas_asum<float_complext>(cublasHandle_t hndl, int n,const float_complext *x, int incx, float *result){
	return cublasScasum(hndl,n,(const cuComplex*) x,incx,result);
}
template<> int cublas_asum<double_complext>(cublasHandle_t hndl, int n,const double_complext *x, int incx, double *result){
	return cublasDzasum(hndl,n,(const cuDoubleComplex*) x,incx,result);
}

//AMIN
template<> int cublas_amin<float>(cublasHandle_t hndl, int n,const float *x, int incx, int *result){
	return cublasIsamin(hndl,n,x,incx,result);
}
template<> int cublas_amin<double>(cublasHandle_t hndl, int n,const double *x, int incx, int *result){
	return cublasIdamin(hndl,n,x,incx,result);
}
template<> int cublas_amin<float_complext>(cublasHandle_t hndl, int n,const float_complext *x, int incx, int *result){
	return cublasIcamin(hndl,n, (const cuComplex* ) x,incx,result);
}
template<> int cublas_amin<double_complext>(cublasHandle_t hndl, int n,const double_complext *x, int incx, int *result){
	return cublasIzamin(hndl,n, (const cuDoubleComplex* ) x,incx,result);
}

//AMAX
template<> int cublas_amax<float>(cublasHandle_t hndl, int n,const float *x, int incx, int *result){
	return cublasIsamax(hndl,n,x,incx,result);
}
template<> int cublas_amax<double>(cublasHandle_t hndl, int n,const double *x, int incx, int *result){
	return cublasIdamax(hndl,n,x,incx,result);
}
template<> int cublas_amax<float_complext>(cublasHandle_t hndl, int n,const float_complext *x, int incx, int *result){
	return cublasIcamax(hndl,n, (const cuComplex* ) x,incx,result);
}
template<> int cublas_amax<double_complext>(cublasHandle_t hndl, int n,const double_complext *x, int incx, int *result){
	return cublasIzamax(hndl,n, (const cuDoubleComplex* ) x,incx,result);
}


template<class T> typename realType<T>::type
nrm2( cuNDArray<T>* arr, int device )
{
	typedef typename realType<T>::type REAL;
	initialize_static_variables();
	REAL ret;

  if (cublas_nrm2<T>( handle[device], arr->get_number_of_elements(),
		    arr->get_data_ptr(), 1,
		   &ret) != CUBLAS_STATUS_SUCCESS )
    {

      throw cuda_error("cuNDANRM2: nrm2 calculation using cublas failed");
    }

  cudaThreadSynchronize();
  return ret;
}

template<class T> T
dot( cuNDArray<T>* arr1, cuNDArray<T>* arr2, int device )
{
	initialize_static_variables();
  T ret;

  if (cublas_dot( handle[device], arr1->get_number_of_elements(),
		   arr1->get_data_ptr(), 1,
		   arr2->get_data_ptr(), 1,
		   &ret) != CUBLAS_STATUS_SUCCESS )
    {
      throw cuda_error("Error CUDA dot");
    }

  cudaThreadSynchronize();
  return ret;
}


template<class T,class R> void
axpy(R a,  cuNDArray<T>* x, cuNDArray<T>* y, int device )
{

	initialize_static_variables();
	T aT = T(a);
  if( cublas_axpy(handle[device], x->get_number_of_elements(),
		  &aT, x->get_data_ptr(), 1,
		  y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS )
    {
	  throw cuda_error("AXPY using cublas failed");

    }

  cudaThreadSynchronize();
}

template<class T> typename realType<T>::type asum(cuNDArray<T>* x,int device){
	initialize_static_variables();
	typename realType<T>::type result;
	if (cublas_asum(handle[device],x->get_number_of_elements(),x->get_data_ptr(),1,&result)){
		throw cuda_error("ASUM using cublas failed");
	}

	return result;
}

template<class T> int amin(cuNDArray<T>* x,int device){
	initialize_static_variables();
	int result;
	if (cublas_amin(handle[device],x->get_number_of_elements(),x->get_data_ptr(),1,&result)){
		throw cuda_error("AMIN using cublas failed");
	}

	return result;
}

template<class T> int amax(cuNDArray<T>* x,int device){
	initialize_static_variables();
	int result;
	if (cublas_amax(handle[device],x->get_number_of_elements(),x->get_data_ptr(),1,&result)){
		throw cuda_error("AMIN using cublas failed");
	}

	return result;
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



