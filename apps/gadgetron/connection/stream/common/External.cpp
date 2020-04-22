#include "External.h"

#include <memory>

#include <boost/asio.hpp>

#include "connection/SocketStreamBuf.h"
#include "log.h"

namespace {

    std::string make_error_message(const std::list<std::string>& errors) {
        std::stringstream error_maker;
        error_maker << "Errors received: " << std::endl;
        for (auto& error : errors) {
            error_maker << error << std::endl;
        }
        return error_maker.str();
    }
}

namespace {
    using namespace Gadgetron::Server::Connection::Stream;

    Remote as_remote(Local, const std::shared_ptr<Configuration> &configuration) {
        return Remote {
                "localhost",
                std::to_string(configuration->context.settings.port)
        };
    }

    Remote as_remote(Remote remote, const std::shared_ptr<Configuration> &) {
        return remote;
    }
}

namespace Gadgetron::Server::Connection::Stream {

    RemoteError::RemoteError(std::list<std::string> messages) : std::runtime_error(make_error_message(messages)), messages(messages) {}

    std::ostream & operator<<(std::ostream &stream, const Local &l) { stream << to_string(l); return stream; }
    std::ostream & operator<<(std::ostream &stream, const Remote &r) { stream << to_string(r); return stream; }

    std::string to_string(const Local &) {
        return "Local";
    }

    std::string to_string(const Remote &remote) {
        return remote.address + ":" + remote.port;
    }

    std::unique_ptr<std::iostream> listen(unsigned short port) {

        boost::asio::io_service service;
        boost::asio::ip::tcp::endpoint peer;
        boost::asio::ip::tcp::endpoint local(boost::asio::ip::tcp::v6(), port);
        boost::asio::ip::tcp::acceptor acceptor(service, local);

        auto socket = std::make_unique<boost::asio::ip::tcp::socket>(service);

        acceptor.accept(*socket, peer);
        acceptor.close();

        GINFO_STREAM("Accepted connection from " << peer.address());

        return Gadgetron::Connection::stream_from_socket(std::move(socket));
    }

    std::unique_ptr<std::iostream> connect(const std::string &address, const std::string &port) {
        return Gadgetron::Connection::remote_stream(address, port);
    }

    std::unique_ptr<std::iostream> connect(const Remote &remote) {
        return connect(remote.address, remote.port);
    }

    std::unique_ptr<std::iostream> connect(const Address &address, std::shared_ptr<Configuration> configuration) {
        return connect(Core::visit([&](auto address) { return as_remote(address, configuration); }, address));
    }
}
