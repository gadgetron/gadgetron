#include "External.h"

#include "connection/Config.h"
#include "connection/SocketStreamBuf.h"

#include "external/Python.h"
#include "external/Matlab.h"

#include <boost/asio.hpp>
#include <boost/algorithm/string.hpp>


using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Stream;

namespace {

    const std::map<std::string, std::function<void(std::string, unsigned short, const Context &)>> modules{
            {"python", start_python_module}
    };

    std::unique_ptr<std::iostream> open_connection(Config::Connect connect, const Context &context) {

        GINFO_STREAM("Connecting to external module on port: " << connect.port);

        return std::unique_ptr<std::iostream>(nullptr);
    }

    std::unique_ptr<std::iostream> open_connection(Config::Execute execute, const Context &context) {

        boost::asio::io_service service;
        boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::v6(), 0);
        boost::asio::ip::tcp::acceptor acceptor(service, endpoint);

        auto port = acceptor.local_endpoint().port();
        boost::algorithm::to_lower(execute.type);

        GINFO_STREAM("Waiting for external module '" << execute.name << "' on port: " << port);

        modules.at(execute.type)(execute.name, port, context);

        auto socket = std::make_unique<boost::asio::ip::tcp::socket>(service);
        acceptor.accept(*socket);
        acceptor.close();

        GINFO_STREAM("Connected to external module '" << execute.name << "' on port: " << port);

        return Gadgetron::Connection::stream_from_socket(std::move(socket));
    }
}

namespace Gadgetron::Server::Connection::Stream {

    External::External(
            const Config::External &config,
            const Core::Context &context,
            Loader &loader
    ) {
        auto connection = boost::apply_visitor([&](auto action) { return open_connection(action, context); }, config.action);
    }

    void External::process(
            Core::InputChannel input,
            Core::OutputChannel output,
            ErrorHandler &error_handler
    ) {

    }

    const std::string &External::name() {
        static const std::string name = "external";
        return name;
    }
}