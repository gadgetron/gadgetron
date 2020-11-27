#pragma once

#include <ostream>
#include <boost/asio/ip/tcp.hpp>
namespace Gadgetron::Server::Info {

    void print_system_information(std::ostream &);


    boost::asio::ip::tcp tcp_protocol();
    std::string ismrmrd_version();

    std::string gadgetron_version();
    std::string gadgetron_build();

    size_t system_memory();
    bool python_support();
    bool matlab_support();

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