#include "External.h"

#include "common/Closer.h"
#include "common/ExternalChannel.h"

#include "connection/SocketStreamBuf.h"
#include "connection/config/Config.h"

#include "external/Python.h"
#include "external/Matlab.h"

#include <boost/asio/use_future.hpp>
#include <boost/algorithm/string.hpp>
#include "system_info.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Nodes;

using tcp = boost::asio::ip::tcp;

namespace {

    const std::map<std::string, std::function<boost::process::child(const Config::Execute &, unsigned short, const Context &)>> modules{
            {"python", start_python_module},
            {"matlab", start_matlab_module}
    };

    void process_input(GenericInputChannel input, std::shared_ptr<ExternalChannel> external) {
        auto closer = make_closer(external);
        for (auto message : input) {
            external->push_message(std::move(message));
        }
    }

    void process_output(OutputChannel output, std::shared_ptr<ExternalChannel> external) {
        while(true) {
            output.push_message(external->pop());
        }
    }
}

namespace Gadgetron::Server::Connection::Nodes {

    void External::monitor_child(
            std::shared_ptr<boost::process::child> child,
            std::shared_ptr<tcp::acceptor> acceptor
    ) {
        child->wait();
        io_service.dispatch([=]() { acceptor->close(); });
    }

    std::shared_ptr<ExternalChannel> External::open_connection(Config::Connect connect, const Context &context) {
        GINFO_STREAM("Connecting to external module on address: " << connect.address << ":" << connect.port);
        return std::make_shared<ExternalChannel>(
                Gadgetron::Connection::remote_stream(connect.address, connect.port),
                serialization,
                configuration
        );
    }

    std::shared_ptr<ExternalChannel> External::open_connection(Config::Execute execute, const Context &context) {

        tcp::endpoint endpoint(Info::tcp_protocol(), 0);
        auto acceptor = std::make_shared<tcp::acceptor>(io_service, endpoint);

        auto port = acceptor->local_endpoint().port();
        boost::algorithm::to_lower(execute.type);

        GINFO_STREAM("Waiting for external module '" << execute.name << "' on port: " << port);

        auto child = std::make_shared<boost::process::child>(modules.at(execute.type)(execute, port, context));

        monitors.child = std::async(
                std::launch::async,
                [=,this](auto child, auto acceptor) { monitor_child(std::move(child), std::move(acceptor)); },
                child,
                acceptor
        );

        auto socket = std::make_unique<tcp::socket>(io_service);
        auto future_socket = acceptor->async_accept(*socket, boost::asio::use_future);

        io_service.run();
        future_socket.get();

        GINFO_STREAM("Connected to external module '" << execute.name << "' on port: " << port);

        auto stream = Gadgetron::Connection::stream_from_socket(std::move(socket));
        auto external_channel = std::make_shared<ExternalChannel>(
                std::move(stream),
                serialization,
                configuration
        );

        return external_channel;
    }

    std::shared_ptr<ExternalChannel> External::open_external_channel(
            const Config::External &config,
            const Context &context
    ) {
        return Core::visit(
                [&, this](auto action) { return this->open_connection(action, context); },
                config.action
        );
    }

    External::External(
            const Config::External &config,
            const Core::StreamContext &context,
            Loader &loader
    ) : serialization(std::make_shared<Serialization>(
                loader.load_default_and_additional_readers(config),
                loader.load_default_and_additional_writers(config)
        )),
        configuration(std::make_shared<Configuration>(
                context,
                config
        )) {
        channel = std::async(
                std::launch::async,
                [=,this](auto config, auto context) { return open_external_channel(config, context); },
                config,
                context
        );
    }

    void External::process(
            InputChannel input,
            OutputChannel output,
            ErrorHandler &error_handler
    ) {
        std::shared_ptr<ExternalChannel> external = channel.get();

        auto input_thread = error_handler.run(
                [=](auto input) { ::process_input(std::move(input), external); },
                std::move(input)
        );

        auto output_thread = error_handler.run(
                [=](auto output) { ::process_output(std::move(output), external); },
                std::move(output)
        );

        input_thread.join(); output_thread.join();
    }

    const std::string &External::name() {
        static const std::string name = "external";
        return name;
    }
}