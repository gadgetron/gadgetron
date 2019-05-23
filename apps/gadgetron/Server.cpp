
#include <boost/asio.hpp>

#include <memory>

#include <Context.h>

#include "log.h"

#include "Server.h"
#include "Connection.h"
#include "connection/SocketStreamBuf.h"

using namespace boost::filesystem;
using namespace Gadgetron::Server;


Server::Server(
        const boost::program_options::variables_map &args
) : args(args) {}

void Server::serve() {

    Gadgetron::Core::Context::Paths paths{args["home"].as<path>(), args["dir"].as<path>()};

    boost::asio::io_service service;
    boost::asio::ip::tcp::endpoint local(boost::asio::ip::tcp::v6(), args["port"].as<unsigned short>());
    boost::asio::ip::tcp::acceptor acceptor(service, local);

    while(true) {
        auto socket = std::make_unique<boost::asio::ip::tcp::socket>(service);
        acceptor.accept(*socket);

        GINFO_STREAM("Accepted connection from: " << socket->remote_endpoint().address());

        Connection::handle(paths, args, Gadgetron::Connection::stream_from_socket(std::move(socket)));
    }
}
