
#include <boost/asio.hpp>

#include <memory>

#include <Context.h>

#include "log.h"

#include "Server.h"
#include "Connection.h"

using namespace boost::filesystem;
using namespace Gadgetron::Server;

Server::Server(boost::asio::io_service &io_service, const boost::program_options::variables_map &args)
    : args_(args) , acceptor_(io_service, tcp::endpoint(tcp::v4(), args["port"].as<unsigned short>())) {
    accept();
}

void Server::accept() {
    acceptor_.async_accept(
            [this](const boost::system::error_code &error, tcp::socket socket) {
                this->connection_handler(error, socket);
                this->accept();
            }
    );
}

void Server::connection_handler(const boost::system::error_code &error, tcp::socket &socket) {

    if (error) {
        GERROR_STREAM("Failed to open connection: " << error);
        return;
    }


    GINFO_STREAM("Accepting connection from: " << socket.remote_endpoint().address());

    auto paths = Gadgetron::Core::Context::Paths(args_["home"].as<path>(), args_["dir"].as<path>());
    Connection::create(paths, socket);
}
