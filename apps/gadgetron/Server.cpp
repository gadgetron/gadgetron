
#include <boost/asio.hpp>

#include <memory>

#include <Context.h>

#include "log.h"

#include "Server.h"
#include "Connection.h"

using namespace boost::filesystem;
using namespace Gadgetron::Server;

#include "connection/SocketStreamBuf.h"

Server::Server(
        boost::asio::io_service &io_service,
        const boost::program_options::variables_map &args
) : args(args),
    acceptor(io_service, tcp::endpoint(tcp::v6(), args["port"].as<unsigned short>())) {
    accept();
}

void Server::accept() {
    // Function accept recurses infinitely: Yes it does.
    socket = std::make_unique<tcp::socket>(acceptor.get_io_service());
    acceptor.async_accept(
            *socket,
            [this](const boost::system::error_code &error) {
                this->connection_handler(error);
                this->accept();
            }
    );
}

void Server::connection_handler(const boost::system::error_code &error) {

    if (error) {
        GERROR_STREAM("Failed to open connection: " << error);
        return;
    }

    GINFO_STREAM("Accepting connection from: " << socket->remote_endpoint().address());

    Gadgetron::Core::Context::Paths paths{args["home"].as<path>(), args["dir"].as<path>()};
    Connection::handle(paths, args, Gadgetron::Connection::stream_from_socket(std::move(socket)));
}
