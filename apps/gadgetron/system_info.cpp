
#include "system_info.h"

#include "gadgetron_config.h"
#include "connection/stream/external/Python.h"
#include "connection/stream/external/Matlab.h"
#include "log.h"


#if defined(_WIN32)
#include <Windows.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))

#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sstream>

#endif
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#if USE_CUDA
// CUDA-C includes
#include <cuda/cuda.h>
#include <cuda/cuda_runtime.h>
#endif

namespace Gadgetron::Server::Info {

    std::string ismrmrd_version() {
        return "0.0.0-alpha";
    }

    std::string gadgetron_version() {
        return GADGETRON_VERSION_STRING;
    }

    std::string gadgetron_build() {
        return GADGETRON_GIT_SHA1_HASH;
    }

    size_t system_memory() {
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
        return (size_t) sysconf(_SC_PHYS_PAGES) *
               (size_t) sysconf(_SC_PAGESIZE);

#endif //Mac

#endif //WIN32
        return 0L;
    }

    bool python_support() {
        return Gadgetron::Server::Connection::Stream::python_available();
    }

    bool matlab_support() {
        return Gadgetron::Server::Connection::Stream::matlab_available();
    }

#if defined USE_CUDA
    namespace CUDA {

        std::string format_version(int version) {
            std::stringstream stream;
            stream << (version / 1000)  << "." << (version % 100) / 10;
            return stream.str();
        }

        int cuda_device_count() {
            int deviceCount = 0;
            auto error = cudaGetDeviceCount(&deviceCount);

            if (error ) return 0;
            return deviceCount;
        }

        bool cuda_support() {
            return 0 < cuda_device_count();
        }

        std::string cuda_driver_version() {
            int version = 0;
            cudaDriverGetVersion(&version);
            return format_version(version);
        }

        std::string cuda_runtime_version() {
            int version = 0;
            cudaRuntimeGetVersion(&version);
            return format_version(version);
        }

        std::string cuda_device_name(int device) {
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, device);

            return properties.name;
        }

        size_t cuda_device_memory(int device) {
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, device);

            return properties.totalGlobalMem;
        }

        std::string cuda_device_capabilities(int device) {
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, device);

            std::stringstream stream;
            stream << properties.major << "." << properties.minor;

            return stream.str();
        }

        void print_cuda_information(std::ostream &os) {
            int device_count = CUDA::cuda_device_count();

            os << "  -- CUDA Support       : YES (" << GADGETRON_CUDA_NVCC_FLAGS << ")" << std::endl;
            os << "    * Number of CUDA capable devices: " << device_count << std::endl;

            for (int dev = 0; dev < device_count; dev++) {
              os << "      - Device " << dev << ": " << CUDA::cuda_device_name(dev) << std::endl;
              os << "         + CUDA Driver Version / Runtime Version: " << CUDA::cuda_driver_version() << "/" << CUDA::cuda_runtime_version() << std::endl;
              os << "         + CUDA Capability Major / Minor version number: " <<  CUDA::cuda_device_capabilities(dev) << std::endl;
              os << "         + Total amount of global GPU memory: " << CUDA::cuda_device_memory(dev) / (1024 * 1024) << " MB" << std::endl;
            }
        }
    }
#else
    namespace CUDA {

        bool cuda_support() {
            return false;
        }

        int cuda_device_count() {
            return 0;
        }

        std::string cuda_driver_version() {
            return "N/A";
        }

        std::string cuda_runtime_version() {
            return "N/A";
        }

        std::string cuda_device_name(int) {
            return "N/A";
        }

        size_t cuda_device_memory(int) {
            return 0;
        }

        std::string cuda_device_capabilities(int) {
            return "N/A";
        }

        void print_cuda_information(std::ostream &os) {
            os << "  -- CUDA Support       : NO" << std::endl;
        }
    }
#endif

    void print_system_information(std::ostream &os) {
        os << "Gadgetron Version Info" << std::endl;
        os << "  -- Version            : " << gadgetron_version().c_str() << std::endl;
        os << "  -- Git SHA1           : " << gadgetron_build().c_str() << std::endl;
        os << "  -- System Memory size : " << system_memory() / (1024 * 1024) << " MB" << std::endl;
        os << "  -- Python Support     : " << (python_support() ? "YES" : "NO") << std::endl;
        os << "  -- Matlab Support     : " << (matlab_support() ? "YES" : "NO") << std::endl;
        CUDA::print_cuda_information(os);
        os << std::endl;
    }
}
