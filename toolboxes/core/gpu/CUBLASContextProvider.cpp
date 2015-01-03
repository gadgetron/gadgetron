#include "CUBLASContextProvider.h"

#include <cuda_runtime_api.h>

#ifdef _WITH_CULA_SUPPORT
#include <cula_lapack_device.h>
#endif


CUBLASContextProvider* CUBLASContextProvider::instance()
{
		if (!instance_) instance_ = new CUBLASContextProvider();
		return instance_;
}

CUBLASContextProvider::~CUBLASContextProvider()
{
	std::map<int, cublasHandle_t>::iterator it = handles_.begin();

	while (it != handles_.end()) {
		if (cudaSetDevice(it->first)!= cudaSuccess) {
		    std::cerr << "Error: unable to set CUDA device." << std::endl;
		}

#ifdef _WITH_CULA_SUPPORT
		culaShutdown();
#endif

		cublasDestroy_v2(it->second);
		it++;
	}

}

cublasHandle_t* CUBLASContextProvider::getCublasHandle(int device_no)
{
	std::map<int, cublasHandle_t>::iterator it;


	//Let's see if we have the handle already:
	it = handles_.find(device_no);

	if (it != handles_.end()) {
		return &handles_[device_no];
	}


	//We don't have the handle yet, let's check if it makes sense to create one

	int number_of_devices = 0;
	if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
	    std::cerr << "Error: unable to query number of CUDA devices.\n" << std::endl;
	    return 0;
	}

	if (number_of_devices == 0) {
	      std::cerr << "Error: No available CUDA devices.\n" << std::endl;
	      return 0;
     }

	  if (device_no >= number_of_devices) {
	      std::cerr << "Requested device number exceeds number of devices." << std::endl;
		  return 0;
	  }

	  //OK, so we are OK to create the handle. Before we do that, let's capture the current cuda device.

	  int current_device_no;
	if (cudaGetDevice(&current_device_no)!= cudaSuccess) {
		 std::cerr << "Error: unable to get current CUDA device.\n" << std::endl;
		      return 0;
	}

	if (current_device_no != device_no) {
		//We must switch context
		if (cudaSetDevice(device_no)!= cudaSuccess) {
		    std::cerr << "Error: unable to set CUDA device." << std::endl;
		      return 0;
		}
	}

	cublasHandle_t handle; // this is a struct pointer

	//GDEBUG_STREAM("*********   CREATING NEW CONTEXT ************" << std::endl);

	if (cublasCreate_v2(&handle) != CUBLAS_STATUS_SUCCESS) {
		std::cerr << "CUBLASContextProvider: unable to create cublas handle\n" << std::endl;
		return 0;
	}

	handles_[device_no] = handle;

#ifdef _WITH_CULA_SUPPORT
	culaStatus s;
	s = culaInitialize();
	if(s != culaNoError) {
		std::cerr << "CUBLASContextProvider: failed to initialize CULA" << std::endl;
		return 0;
	}
#endif

	if (current_device_no != device_no) {
		//We must switch context back
		if (cudaSetDevice(current_device_no)!= cudaSuccess) {
		   std::cerr << "Error: unable to set CUDA device.\n" << std::endl;
		    return 0;
		}
	}

	return &handles_[device_no];
}


CUBLASContextProvider* CUBLASContextProvider::instance_ = 0;

