
#include <boost/algorithm/string/join.hpp>
#include "Handlers.h"

#include "system_info.h"

#include "io/primitives.h"
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

    void initialize_with_default_queries(std::map<std::string, std::function<std::string()>> &answers) {
        answers["ismrmrd::version"]              = Info::ismrmrd_version;
        answers["gadgetron::version"]            = Info::gadgetron_version;
        answers["gadgetron::build"]              = Info::gadgetron_build;
        answers["gadgetron::info"]               = gadgetron_info;
        answers["gadgetron::info::memory"]       = []() { return std::to_string(Info::system_memory()); };
        answers["gadgetron::info::python"]       = []() { return std::to_string(Info::python_support()); };
        answers["gadgetron::info::matlab"]       = []() { return std::to_string(Info::matlab_support()); };
        answers["gadgetron::info::cuda"]         = []() { return std::to_string(Info::CUDA::cuda_support()); };
        answers["gadgetron::cuda::devices"]      = []() { return std::to_string(Info::CUDA::cuda_device_count()); };
        answers["gadgetron::cuda::driver"]       = Info::CUDA::cuda_driver_version;
        answers["gadgetron::cuda::runtime"]      = Info::CUDA::cuda_runtime_version;
        answers["gadgetron::cuda::memory"]       = cuda_memory;
        answers["gadgetron::cuda::capabilities"] = cuda_capabilities;
    }
}

namespace Gadgetron::Server::Connection::Handlers {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::IO;

    QueryHandler::QueryHandler() {
        initialize_with_default_queries(answers);
    }

    void QueryHandler::handle(std::istream &stream, Gadgetron::Core::OutputChannel& channel) {

        auto reserved = read<uint64_t>(stream);
        auto corr_id  = read<uint64_t>(stream);
        auto query    = read_string_from_stream<uint64_t>(stream);

        if (reserved) {
            throw std::runtime_error("Unsupported value in reserved bytes.");
        }

        channel.push(Response(corr_id, answers.at(query)()));
    }


    ErrorProducingHandler::ErrorProducingHandler(std::string message)
    : message(std::move(message)) {}


    void ErrorProducingHandler::handle(std::istream &, Gadgetron::Core::OutputChannel&) {
        throw std::runtime_error(message);
    }

    CloseHandler::CloseHandler(std::function<void()> callback) : callback(std::move(callback)) {}

    void CloseHandler::handle(std::istream &stream, Gadgetron::Core::OutputChannel& ) { callback(); }
}

