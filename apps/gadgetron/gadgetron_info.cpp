//#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
//#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"
#include "ace/DLL_Manager.h"
//#include "ace/OS_NS_netdb.h"

#include "gadgetron_config.h"
#include "Gadget.h"

#include <iostream>


using namespace Gadgetron;

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

size_t get_system_memory_size()
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

int main(int argc, char** argv)
{
  GDEBUG_STREAM("Gadgetron Version Info" << std::endl);
  GDEBUG_STREAM("  -- Version            : " << GADGETRON_VERSION_STRING << std::endl);
  GDEBUG_STREAM("  -- Git SHA1           : " << GADGETRON_GIT_SHA1_HASH << std::endl);
  GDEBUG_STREAM("  -- System Memory size : " << get_system_memory_size()/(1024*1024) << " MB" << std::endl);

#if defined COMPILING_WITH_PYTHON_SUPPORT
  GDEBUG_STREAM("  -- Python Support     : YES" << std::endl);
#else
  GDEBUG_STREAM("  -- Python Support     : NO" << std::endl); 
#endif

#if defined USE_CUDA
  GDEBUG_STREAM("  -- CUDA Support       : YES (" << GADGETRON_CUDA_NVCC_FLAGS << ")" << std::endl);
  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
  
  if (error_id != cudaSuccess) {
    GDEBUG_STREAM("    * Unable to get device count" << std::endl);
  } else {
    GDEBUG_STREAM("    * Number of CUDA capable devices: " << deviceCount << std::endl);
    if (deviceCount) {
      
      int dev, driverVersion = 0, runtimeVersion = 0;
      for (dev = 0; dev < deviceCount; ++dev) {
	
        cudaSetDevice(dev);
	
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&runtimeVersion);

	GDEBUG_STREAM("      - Device " << dev << ": " << deviceProp.name << std::endl); 
	std::cout << "         + CUDA Driver Version / Runtime Version: " 
		  <<  (driverVersion/1000)  << "." << (driverVersion%100)/10 << "/" 
		  <<  (runtimeVersion/1000) << "." << (runtimeVersion%100)/10 << std::endl;
	GDEBUG_STREAM("         + CUDA Capability Major/Minor version number: " <<  deviceProp.major << "." << deviceProp.minor << std::endl);
	GDEBUG_STREAM("         + Total amount of global GPU memory: " << (float)deviceProp.totalGlobalMem/1048576.0f << " MB" << std::endl);
      }
    }
  }
#else
  GDEBUG_STREAM("  -- CUDA Support       : NO" << std::endl); 
#endif
    
  GDEBUG_STREAM(std::endl);

  if (argc == 1) {
    return 0;
  }

  if ((argc == 2) || (argc > 3)) {
    GDEBUG_STREAM("Invalid number of arguments (" << argc -1 << ")." << std::endl);
    GDEBUG_STREAM("Usage (gadget library info):  " << argc << std::endl);
    GDEBUG_STREAM(" -- gadgetron_info <SHARED LIB> <GADGET_INFO>" << std::endl);
    return -1; 
  }

  const char* DLL = argv[1];
  const char* component_name = argv[2];

  //We must be investigating a certain gadget
  GDEBUG_STREAM("Examining Gadget (SHARED LIB): " << component_name << " (" << DLL << ")" << std::endl);

  //Attempt to load Gadget
  //ACE_DLL_Manager* dllmgr = ACE_DLL_Manager::instance();
  
  ACE_DLL_Handle dll;// = 0;
  ACE_SHLIB_HANDLE dll_handle = 0;
  
  ACE_TCHAR dllname[1024];
#if defined(WIN32) && defined(_DEBUG)
  ACE_OS::sprintf(dllname, "%s%sd",ACE_DLL_PREFIX, DLL);
#else
  ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);
#endif

  ACE_TCHAR factoryname[1024];
  ACE_OS::sprintf(factoryname, "make_%s", component_name);
  
  if (dll.open(dllname, ACE_DEFAULT_SHLIB_MODE, dll_handle )) {
    GDEBUG_STREAM("Failed to load DLL (" << DLL << "), Possible reasons:" << std::endl);
    GDEBUG_STREAM("   - Name of DLL is wrong" << std::endl);
    GDEBUG_STREAM("   - Path of DLL is not in your DLL search path (LD_LIBRARY_PATH on Unix)" << std::endl);
    GDEBUG_STREAM("   - Path of other DLLs that this DLL depends on is not in the search path" << std::endl);
    GDEBUG_STREAM("" << std::endl);
    GDEBUG_STREAM("Set environment variable ACE_DEBUG=1 to get more information" << std::endl << std::endl); 
    return 0;
  } 

  //Function pointer
  typedef Gadget* (*ComponentCreator) (void);

  void *void_ptr = dll.symbol (factoryname);
  ptrdiff_t tmp = reinterpret_cast<ptrdiff_t> (void_ptr);
  ComponentCreator cc = reinterpret_cast<ComponentCreator> (tmp);
  
  if (cc == 0) {
    GDEBUG_STREAM("Failed to load factory (" << factoryname << ") from DLL (" << dllname << ")" << std::endl);
    return -1;
  }
  
  Gadget* g = cc();
  if (!g) {
    GDEBUG_STREAM("Failed to create component using factory" << std::endl);
    return 0;
  }

  GDEBUG_STREAM("  -- Gadget compiled against Gadgetron version " << g->get_gadgetron_version() << std::endl);

  delete g;

  return 0;
}
