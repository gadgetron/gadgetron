#ifndef GADGETRONSYSTEMINFO_H
#define GADGETRONSYSTEMINFO_H

#include "gadgetron_config.h"

#include <ostream>

#if defined(_WIN32)
#include <Windows.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#endif
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#if USE_CUDA
// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime.h>
#endif


namespace Gadgetron
{
  
  inline size_t get_system_memory_size()
  {
#if defined(_WIN32)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx( &status );
    return (size_t)status.ullTotalPhys;
#else //Unix variant

#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64)) //Mac
    int mib[2];
    mib[0] = CTL_HW;

#if defined(HW_MEMSIZE)
    mib[1] = HW_MEMSIZE;
#elif defined(HW_PHYSMEM64)
    mib[1] = HW_PHYSMEM64;
#endif

    int64_t size = 0;
    size_t len = sizeof( size );
    if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
      return (size_t)size;
    return 0L;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE) //Linux
    return (size_t)sysconf( _SC_PHYS_PAGES ) *
      (size_t)sysconf( _SC_PAGESIZE );
  
#endif //Mac

#endif //WIN32
    return 0L;
  }

  

  inline void print_system_information(std::ostream& os)
  {

    os << "Gadgetron Version Info" << std::endl;
    os << "  -- Version            : " << GADGETRON_VERSION_STRING << std::endl;
    os << "  -- Git SHA1           : " << GADGETRON_GIT_SHA1_HASH << std::endl;
    os << "  -- System Memory size : " << get_system_memory_size()/(1024*1024) << " MB" << std::endl;

#if defined COMPILING_WITH_PYTHON_SUPPORT
    os << "  -- Python Support     : YES" << std::endl;
#else
    os << "  -- Python Support     : NO" << std::endl; 
#endif

#if defined USE_MATLAB
    os << "  -- Matlab Support     : YES" << std::endl;
#else
    os << "  -- Matlab Support     : NO" << std::endl;
#endif


#if defined USE_CUDA
    os << "  -- CUDA Support       : YES (" << GADGETRON_CUDA_NVCC_FLAGS << ")" << std::endl;
    int deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
  
    if (error_id != cudaSuccess) {
      os << "    * Unable to get device count" << std::endl;
    } else {
      os << "    * Number of CUDA capable devices: " << deviceCount << std::endl;
      if (deviceCount) {
      
	int dev, driverVersion = 0, runtimeVersion = 0;
	for (dev = 0; dev < deviceCount; ++dev) {
	
	  cudaSetDevice(dev);
	
	  cudaDeviceProp deviceProp;
	  cudaGetDeviceProperties(&deviceProp, dev);
	  cudaDriverGetVersion(&driverVersion);
	  cudaRuntimeGetVersion(&runtimeVersion);

	  os << "      - Device " << dev << ": " << deviceProp.name << std::endl; 
	  os << "         + CUDA Driver Version / Runtime Version: " 
	     <<  (driverVersion/1000)  << "." << (driverVersion%100)/10 << "/" 
	     <<  (runtimeVersion/1000) << "." << (runtimeVersion%100)/10 << std::endl;
	  os << "         + CUDA Capability Major/Minor version number: " <<  deviceProp.major << "." << deviceProp.minor << std::endl;
	  os << "         + Total amount of global GPU memory: " << (float)deviceProp.totalGlobalMem/1048576.0f << " MB" << std::endl;
	}
      }
    }
#else
    os << "  -- CUDA Support       : NO" << std::endl; 
#endif
  
    os << std::endl;
  
  }
}

#endif //GADGETRONSYSTEMINFO_H
