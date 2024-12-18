#pragma once

#include <ostream>
#include <boost/filesystem/path.hpp>

namespace Gadgetron::Main::Info {

    void print_system_information(std::ostream &);

    std::string pingvin_version();
    std::string pingvin_build();

    size_t system_memory();
    bool python_support();
    bool matlab_support();

    const boost::filesystem::path default_pingvin_home();

    namespace CUDA {
        bool cuda_support();
        int cuda_device_count();
        std::string cuda_driver_version();
        std::string cuda_runtime_version();
        std::string cuda_device_name(int);
        size_t cuda_device_memory(int);
        std::string cuda_device_capabilities(int);
    }
}