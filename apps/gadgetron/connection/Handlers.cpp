
#include <boost/algorithm/string/join.hpp>
#include "Handlers.h"

#include "gadgetron_system_info.h"

#include "readers/Primitives.h"
#include "Response.h"

namespace {

    using namespace Gadgetron::Server;

    std::string gadgetron_info() {
        std::stringstream stream;
        Info::print_system_information(stream);
        return stream.str();
    }

    std::vector<std::string> for_each_device(std::function<std::string(int)> fn) {
        std::vector<std::string> results{};
        for (int device = 0; device < Info::CUDA::cuda_device_count(); device++) {
            results.push_back(fn(device));
        }
        return results;
    }

    std::string cuda_capabilities() {
        return boost::algorithm::join(for_each_device(Info::CUDA::cuda_device_capabilities), ";");
    }

    std::string cuda_memory() {
        return boost::algorithm::join(
                for_each_device(
                        [](int device) {
                            return std::to_string(Info::CUDA::cuda_device_memory(device));
                        }
                ),
                ";"
        );
    }

    void initialize_with_default_queries(std::map<std::string, std::function<std::string()>> &handlers) {

        handlers["ismrmrd::version"]              = Info::ismrmrd_version;
        handlers["gadgetron::version"]            = Info::gadgetron_version;
        handlers["gadgetron::build"]              = Info::gadgetron_build;
        handlers["gadgetron::info"]               = gadgetron_info;
        handlers["gadgetron::info::memory"]       = []() { return std::to_string(Info::system_memory()); };
        handlers["gadgetron::info::python"]       = []() { return std::to_string(Info::python_support()); };
        handlers["gadgetron::info::matlab"]       = []() { return std::to_string(Info::matlab_support()); };
        handlers["gadgetron::info::cuda"]         = []() { return std::to_string(Info::CUDA::cuda_support()); };
        handlers["gadgetron::cuda::devices"]      = []() { return std::to_string(Info::CUDA::cuda_device_count()); };
        handlers["gadgetron::cuda::driver"]       = Info::CUDA::cuda_driver_version;
        handlers["gadgetron::cuda::runtime"]      = Info::CUDA::cuda_runtime_version;
        handlers["gadgetron::cuda::memory"]       = cuda_memory;
        handlers["gadgetron::cuda::capabilities"] = cuda_capabilities;
    }
}

namespace Gadgetron::Server::Connection::Handlers {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;

    QueryHandler::QueryHandler(std::shared_ptr<Gadgetron::Core::OutputChannel> channel)
    : channel(std::move(channel)) {
        initialize_with_default_queries(handlers);
    }

    void QueryHandler::handle(std::istream &stream) {

        auto reserved = read_t<uint64_t>(stream);
        auto corr_id  = read_t<uint64_t>(stream);
        auto query    = read_string_from_stream<uint64_t>(stream);

        if (reserved) {
            throw std::runtime_error("Unsupported value in reserved bytes.");
        }

        std::string response = handlers.at(query)();

        channel->push(std::make_unique<Response>(corr_id, response));
    }


    ErrorProducingHandler::ErrorProducingHandler(std::string message)
    : message(std::move(message)) {}


    void ErrorProducingHandler::handle(std::istream &) {
        throw std::runtime_error(message);
    }
}

