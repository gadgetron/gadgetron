#include "cudaDeviceManager.h"
#include "check_CUDA.h"
#include "cuNDArray_blas.h"

#include <boost/thread/mutex.hpp>
#include <boost/shared_array.hpp>
#include <cuda_runtime_api.h>
#include <stdlib.h>
#include <sstream>

 std::string Gadgetron::gadgetron_getCusparseErrorString(cusparseStatus_t err)
  {
    switch (err)
    {
    case CUSPARSE_STATUS_NOT_INITIALIZED:
      return "NOT INITIALIZED";
    case CUSPARSE_STATUS_ALLOC_FAILED:
      return "ALLOC FAILED";
    case CUSPARSE_STATUS_INVALID_VALUE:
      return "INVALID VALUE";
    case CUSPARSE_STATUS_ARCH_MISMATCH:
      return "ARCH MISMATCH";
    case CUSPARSE_STATUS_MAPPING_ERROR:
      return "MAPPING ERROR";
    case CUSPARSE_STATUS_EXECUTION_FAILED:
      return "EXECUTION FAILED";
    case CUSPARSE_STATUS_INTERNAL_ERROR:
      return "INTERNAL ERROR";
    case CUSPARSE_STATUS_SUCCESS:
      return "SUCCES";
    case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      return "MATRIX TYPE NOT SUPPORTED";
    default:
      return "UNKNOWN CUSPARSE ERROR";
    }
  }
namespace Gadgetron
{

 

  static boost::shared_array<boost::mutex> _mutex;
  static boost::shared_array<boost::mutex> _sparseMutex;

  cudaDeviceManager *cudaDeviceManager::_instance = 0;

  cudaDeviceManager::cudaDeviceManager()
  {

    // This constructor is executed only once for a singleton
    //

    atexit(&CleanUp);

    if (auto res = cudaGetDeviceCount(&_num_devices); res != cudaSuccess)
    {
      _num_devices = 0;
      throw cuda_error(res);
    }

    _mutex = boost::shared_array<boost::mutex>(new boost::mutex[_num_devices]);
    _sparseMutex = boost::shared_array<boost::mutex>(new boost::mutex[_num_devices]);

    int old_device;
    if (auto res = cudaGetDevice(&old_device); res != cudaSuccess)
    {
      throw cuda_error(res);
    }

    _total_global_mem = std::vector<size_t>(_num_devices, 0);
    _shared_mem_per_block = std::vector<size_t>(_num_devices, 0);
    _warp_size = std::vector<int>(_num_devices, 0);
    _max_blockdim = std::vector<int>(_num_devices, 0);
    _max_griddim = std::vector<int>(_num_devices, 0);
    _major = std::vector<int>(_num_devices, 0);
    _minor = std::vector<int>(_num_devices, 0);
    _handle = std::vector<cublasHandle_t>(_num_devices, (cublasContext *)0x0);
    _sparse_handle = std::vector<cusparseHandle_t>(_num_devices, (cusparseHandle_t)0x0);

    for (int device = 0; device < _num_devices; device++)
    {

      if (auto res = cudaSetDevice(device); res != cudaSuccess)
      {
        throw cuda_error(res);
      }

      cudaDeviceProp deviceProp;

      if (auto res = cudaGetDeviceProperties(&deviceProp, device); res != cudaSuccess)
      {
        throw cuda_error(res);
      }

      _total_global_mem[device] = deviceProp.totalGlobalMem;
      _shared_mem_per_block[device] = deviceProp.sharedMemPerBlock;
      _warp_size[device] = deviceProp.warpSize;
      _max_blockdim[device] = deviceProp.maxThreadsDim[0];
      _max_griddim[device] = deviceProp.maxGridSize[0];
      _major[device] = deviceProp.major;
      _minor[device] = deviceProp.minor;
    }

    if (auto res = cudaSetDevice(old_device); res != cudaSuccess)
    {
      throw cuda_error(res);
    }
  }

  cudaDeviceManager::~cudaDeviceManager()
  {

    for (int device = 0; device < _num_devices; device++)
    {
      if (_handle[device] != NULL)
        cublasDestroy(_handle[device]);
      if (_sparse_handle[device] != NULL)
        cusparseDestroy(_sparse_handle[device]);
    }
  }

  size_t cudaDeviceManager::total_global_mem()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return _total_global_mem[device];
  }

  size_t cudaDeviceManager::shared_mem_per_block()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return _shared_mem_per_block[device];
  }

  int cudaDeviceManager::max_blockdim()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return _max_blockdim[device];
  }

  int cudaDeviceManager::max_griddim()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return _max_griddim[device];
  }

  int cudaDeviceManager::warp_size()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return _warp_size[device];
  }

  int cudaDeviceManager::major_version()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return _major[device];
  }

  int cudaDeviceManager::minor_version()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return _minor[device];
  }

  size_t cudaDeviceManager::getFreeMemory()
  {
    size_t free, total;
    CUDA_CALL(cudaMemGetInfo(&free, &total));
    return free;
  }

  size_t cudaDeviceManager::getTotalMemory()
  {
    size_t free, total;
    CUDA_CALL(cudaMemGetInfo(&free, &total));
    return total;
  }

  size_t cudaDeviceManager::getFreeMemory(int device)
  {
    int oldDevice;
    CUDA_CALL(cudaGetDevice(&oldDevice));
    CUDA_CALL(cudaSetDevice(device));
    size_t ret = getFreeMemory();
    CUDA_CALL(cudaSetDevice(oldDevice));
    return ret;
  }

  size_t cudaDeviceManager::getTotalMemory(int device)
  {
    int oldDevice;
    CUDA_CALL(cudaGetDevice(&oldDevice));
    CUDA_CALL(cudaSetDevice(device));
    size_t ret = getTotalMemory();
    CUDA_CALL(cudaSetDevice(oldDevice));
    return ret;
  }

  cudaDeviceManager *cudaDeviceManager::Instance()
  {
    if (_instance == 0)
      _instance = new cudaDeviceManager;
    return _instance;
  }

  cublasHandle_t cudaDeviceManager::lockHandle()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return lockHandle(device);
  }

  cublasHandle_t cudaDeviceManager::lockHandle(int device)
  {
    _mutex[device].lock();
    if (_handle[device] == NULL)
    {
      cublasStatus_t ret = cublasCreate(&_handle[device]);
      if (ret != CUBLAS_STATUS_SUCCESS)
      {
        std::stringstream ss;
        ss << "Error: unable to create cublas handle for device " << device << " : ";
        ss << gadgetron_getCublasErrorString(ret) << std::endl;
        throw cuda_error(ss.str());
      }
      cublasSetPointerMode(_handle[device], CUBLAS_POINTER_MODE_HOST);
    }
    return _handle[device];
  }

  void cudaDeviceManager::unlockHandle()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return unlockHandle(device);
  }

  void cudaDeviceManager::unlockHandle(int device)
  {
    _mutex[device].unlock();
  }

  cusparseHandle_t cudaDeviceManager::lockSparseHandle()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return lockSparseHandle(device);
  }

  cusparseHandle_t cudaDeviceManager::lockSparseHandle(int device)
  {
    _sparseMutex[device].lock();
    if (_sparse_handle[device] == NULL)
    {
      cusparseStatus_t ret = cusparseCreate(&_sparse_handle[device]);
      if (ret != CUSPARSE_STATUS_SUCCESS)
      {
        std::stringstream ss;
        ss << "Error: unable to create cusparse handle for device " << device << " : ";
        ss << gadgetron_getCusparseErrorString(ret) << std::endl;
        throw cuda_error(ss.str());
      }
      cusparseSetPointerMode(_sparse_handle[device], CUSPARSE_POINTER_MODE_HOST);
      //cublasSetPointerMode( _handle[device], CUBLAS_POINTER_MODE_HOST );
    }
    return _sparse_handle[device];
  }

  void cudaDeviceManager::unlockSparseHandle()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return unlockSparseHandle(device);
  }

  void cudaDeviceManager::unlockSparseHandle(int device)
  {
    _sparseMutex[device].unlock();
  }

  int cudaDeviceManager::getCurrentDevice()
  {
    int device;
    CUDA_CALL(cudaGetDevice(&device));
    return device;
  }

  int cudaDeviceManager::getTotalNumberOfDevice()
  {
    int number_of_devices;
    CUDA_CALL(cudaGetDeviceCount(&number_of_devices));
    return number_of_devices;
  }

  void cudaDeviceManager::CleanUp()
  {
    delete _instance;
    _instance = 0;
  }
} // namespace Gadgetron
