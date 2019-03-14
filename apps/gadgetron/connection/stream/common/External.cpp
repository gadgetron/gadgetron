#include "External.h"

#include <memory>

#include <boost/asio.hpp>

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

    void set_exception_bits(boost::asio::ip::tcp::iostream &stream) {
        stream.exceptions(
                boost::asio::ip::tcp::iostream::failbit |
                boost::asio::ip::tcp::iostream::badbit |
                boost::asio::ip::tcp::iostream::eofbit
        );
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

        boost::asio::io_context context;
        boost::asio::ip::tcp::endpoint peer;
        boost::asio::ip::tcp::endpoint local(boost::asio::ip::tcp::v6(), port);
        boost::asio::ip::tcp::acceptor acceptor(context, local);

        auto stream = std::make_unique<boost::asio::ip::tcp::iostream>();

        acceptor.accept(*stream->rdbuf(), peer);
        acceptor.close();

        GINFO_STREAM("Accepted connection from " << peer.address());

        set_exception_bits(*stream);

        return std::move(stream);
    }

    std::unique_ptr<std::iostream> connect(std::string address, std::string service) {

        auto stream = std::make_unique<boost::asio::ip::tcp::iostream>(boost::asio::ip::tcp::v6(), address, service);

        GINFO_STREAM("Connected to " << address << ":" << service);

        set_exception_bits(*stream);

        return std::move(stream);
    }

    std::unique_ptr<std::iostream> connect(const Remote &remote) {
        return connect(remote.address, remote.port);
    }
}
