
#include "system_info.h"

#include "pingvin_config.h"
#include "log.h"
#include "Process.h"

#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#include <pwd.h>
#include <sstream>
#include <memory>
#include <boost/filesystem.hpp>

#if defined(BSD)
#include <sys/sysctl.h>
#endif

#if USE_CUDA
// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime.h>
#endif


namespace Gadgetron::Main::Info {

    std::string pingvin_version() {
        return PINGVIN_VERSION_STRING;
    }

    std::string pingvin_build() {
        return PINGVIN_GIT_SHA1_HASH;
    }

    size_t system_memory() {

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

        return 0L;
    }

    namespace CUDA {
#if defined USE_CUDA
        std::string format_version(int version) {
            std::stringstream stream;
            stream << (version / 1000)  << "." << (version % 100) / 10;
            return stream.str();
        }

        int cuda_device_count() {
            int deviceCount = 0;
            auto error = cudaGetDeviceCount(&deviceCount);
            GINFO("CUDA DEVICE COUNT %d and error number %u\n", deviceCount, error);
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

            os << "  -- CUDA Support       : YES" << std::endl;
            os << "  -- NVCC Flags         : " << PINGVIN_CUDA_NVCC_FLAGS << std::endl;
            os << "  -- CUDA Device count  : " << device_count << std::endl;

            for (int dev = 0; dev < device_count; dev++) {
              os << "      - Device " << dev << ": " << CUDA::cuda_device_name(dev) << std::endl;
              os << "         + CUDA Driver Version / Runtime Version: " << CUDA::cuda_driver_version() << "/" << CUDA::cuda_runtime_version() << std::endl;
              os << "         + CUDA Capability Major / Minor version number: " <<  CUDA::cuda_device_capabilities(dev) << std::endl;
              os << "         + CUDA Device Memory size: " << std::to_string(CUDA::cuda_device_memory(dev) / (1024 * 1024)) << " MB" << std::endl;
            }
        }
#else
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
#endif
    } // namespace CUDA

    void print_system_information(std::ostream &os) {
        os << "Pingvin Version Info" << std::endl;
        os << "  -- Version            : " << pingvin_version().c_str() << std::endl;
        os << "  -- Git SHA1           : " << pingvin_build().c_str() << std::endl;
        os << "  -- System Memory size : " << std::to_string(system_memory() / (1024 * 1024)) << " MB" << std::endl;
        CUDA::print_cuda_information(os);
        os << std::endl;
    }

    namespace {

    #if defined __APPLE__
        #include <mach-o/dyld.h>

        std::string get_executable_path() {
            char path[PATH_MAX];
            char resolved[PATH_MAX];
            uint32_t size = sizeof(path);

            if ((_NSGetExecutablePath(path, &size) == 0) && (realpath(path, resolved) != NULL)) {
                return std::string(resolved);
            } else {
                throw std::runtime_error("Could not determine location of Pingvin binary.");
            }
        }
    #else
        std::string get_executable_path(size_t buffer_size = 1024) {

            auto buffer = std::make_unique<char[]>(buffer_size);

            ssize_t len = readlink("/proc/self/exe", buffer.get(), buffer_size);

            if (len < 0) {
                throw std::runtime_error("Failed to read /proc/self/exe - cannot determine Pingvin binary path.");
            }

            if (size_t(len) == buffer_size) {
                // Allocated buffer was probably too small. Try again with a bigger buffer.
                return get_executable_path(buffer_size * 2);
            }

            return {buffer.get(), size_t(len)};
        }
    #endif
    } // namespace

    const boost::filesystem::path default_pingvin_home() {

        const char *home = std::getenv("PINGVIN_HOME");

        if (home != nullptr) {
            return boost::filesystem::path(home);
        }

        boost::filesystem::path executable_path = get_executable_path();

        GDEBUG_STREAM("Executable path: " << executable_path);

        boost::filesystem::path pingvin_home = executable_path
                .parent_path()
                .parent_path();

        GDEBUG_STREAM("Default Pingvin home: " << pingvin_home);

        return pingvin_home;
    }

} // namespace
