#pragma once

#include "gpucore_export.h"

#include <vector>
#include <cublas_v2.h>
#include <cusparse.h>

#include <library_types.h>
#include "complext.h"

namespace Gadgetron{

  class EXPORTGPUCORE cudaDeviceManager 
  {
  public:

    // This class is used as a singleton.
    // Use Instance() to access the public member functions.
    //

    static cudaDeviceManager* Instance();
    
    // Public member functions.
    // If the function does not take a device id, it will use the current device.
    //

    inline size_t total_global_mem(int device){ return _total_global_mem[device]; }
    inline size_t shared_mem_per_block(int device){ return _shared_mem_per_block[device]; }
    inline int warp_size(int device){ return _warp_size[device]; }
    inline int max_blockdim(int device){ return _max_blockdim[device]; }
    inline int max_griddim(int device){ return _max_griddim[device]; }
    inline int major_version(int device){ return _major[device]; }
    inline int minor_version(int device){ return _minor[device]; }

    size_t total_global_mem();
    size_t shared_mem_per_block();
    int major_version();
    int minor_version();
    int warp_size();
    int max_blockdim();
    int max_griddim();

    int getCurrentDevice();

    int getTotalNumberOfDevice();

    size_t getFreeMemory();
    size_t getFreeMemory(int device);

    size_t getTotalMemory();
    size_t getTotalMemory(int device);

    // Access to Cublas is protected by a mutex
    // Despite what the Cublas manual claims, we have not found it thread safe.

    cublasHandle_t lockHandle();
    cublasHandle_t lockHandle(int device);

    void unlockHandle();
    void unlockHandle(int device);

    cusparseHandle_t lockSparseHandle();
    cusparseHandle_t lockSparseHandle(int device);

    void unlockSparseHandle();
    void unlockSparseHandle(int device);


  private:

    // Use the Instance() method to access the singleton
    //

    cudaDeviceManager();
    ~cudaDeviceManager();

    static void CleanUp();
    
    int _num_devices;
    std::vector<size_t> _total_global_mem; // in bytes
    std::vector<size_t> _shared_mem_per_block; // in bytes
    std::vector<int> _warp_size;
    std::vector<int> _max_blockdim;
    std::vector<int> _max_griddim;
    std::vector<int> _major;
    std::vector<int> _minor;
    std::vector<cublasHandle_t> _handle;
    std::vector<cusparseHandle_t> _sparse_handle;
    static cudaDeviceManager * _instance;
  };

  template<class T>
  struct cudaDataType {
  };

  template<>
  struct cudaDataType<float>{
    static constexpr cudaDataType_t value = CUDA_R_32F;
  };

  template<>
  struct cudaDataType<double>{
    static constexpr cudaDataType_t value = CUDA_R_64F;
  };

  template<>
  struct cudaDataType<complext<float>>{
    static constexpr cudaDataType_t value = CUDA_C_32F;
  };

  template<>
  struct cudaDataType<complext<double>>{
    static constexpr cudaDataType_t value = CUDA_C_64F;
  };


  template<class T>
  constexpr cudaDataType_t cuda_datatype(){return cudaDataType<T>::value;}

	std::string gadgetron_getCusparseErrorString(cusparseStatus_t err);



}



