
#include "Query.h"

#include "gadgetron_system_info.h"

namespace {

    using namespace Gadgetron::Server;

    std::string concat_for_cuda_devices(const std::function<std::string(int)> &fn, int device_index = 0) {
        if (device_index <= Info::CUDA::cuda_device_count()) {
            return fn(device_index);
        }
        return fn(device_index) + ";" + concat_for_cuda_devices(fn, device_index + 1);
    }

}

namespace Gadgetron::Server::Query {

    GadgetronHandler::GadgetronHandler() : Handler() {

        handlers["gadgetron::version"]       = Info::gadgetron_version;
        handlers["gadgetron::build"]         = Info::gadgetron_build;

        handlers["gadgetron::info"] = []() {
            std::stringstream stream;
            Info::print_system_information(stream);
            return stream.str();
        };

        handlers["gadgetron::info::memory"]  = []() { return std::to_string(Info::system_memory()); };
        handlers["gadgetron::info::python"]  = []() { return std::to_string(Info::python_support()); };
        handlers["gadgetron::info::matlab"]  = []() { return std::to_string(Info::matlab_support()); };
        handlers["gadgetron::info::cuda"]    = []() { return std::to_string(Info::CUDA::cuda_support()); };

        handlers["gadgetron::cuda::driver"]  = Info::CUDA::cuda_driver_version;
        handlers["gadgetron::cuda::runtime"] = Info::CUDA::cuda_runtime_version;

        handlers["gadgetron::cuda::capabilities"] = []() {
            return concat_for_cuda_devices(Info::CUDA::cuda_device_capabilities);
        };
        handlers["gadgetron::cuda::memory"] = []() {
            return concat_for_cuda_devices(
                    [](int device) { return std::to_string(Info::CUDA::cuda_device_memory(device)); }
            );
        };
    }

    bool GadgetronHandler::accepts(const std::string &query) {
        return handlers.count(query) != 0;
    }

    std::string GadgetronHandler::handle(const std::string &query) {
        return handlers.at(query)();
    }

    bool ISMRMRDHandler::accepts(const std::string &query) {
        return query == "ismrmrd::version";
    }

    std::string ISMRMRDHandler::handle(const std::string &query) {
        return Info::ismrmrd_version();
    }

    ConfigHandler::ConfigHandler(std::future<std::stringstream> &&raw_config_future)
    : config(std::move(raw_config_future)) {}

    bool ConfigHandler::accepts(const std::string &query) {
        return query == "gadgetron::config";
    }

    std::string ConfigHandler::handle(const std::string &query) {
        return config.get().str();
    }
}