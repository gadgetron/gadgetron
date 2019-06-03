
#include <boost/asio.hpp>

#include <memory>

#include <Context.h>

#include "log.h"

#include "Server.h"
#include "Connection.h"

using namespace boost::filesystem;
using namespace Gadgetron::Server;
#include "connection/SocketStreamBuf.h"

Server::Server(boost::asio::io_service &io_service, const boost::program_options::variables_map &args)
    : args_(args) , acceptor_(io_service, tcp::endpoint(tcp::v6(), args["port"].as<unsigned short>())) {
    accept();
}

void Server::accept() {
    // Function accept recurses infinitely: Yes it does.
#if(BOOST_VERSION >= 107000)
    socket = std::make_unique<tcp::socket>(acceptor_.get_executor());
#else
    socket = std::make_unique<tcp::socket>(acceptor_.get_io_service());
#endif
    acceptor_.async_accept(
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

    auto paths = Gadgetron::Core::Context::Paths{args_["home"].as<path>(), args_["dir"].as<path>()};
    Connection::handle(paths,Gadgetron::Connection::stream_from_socket(std::move(socket)));
}
