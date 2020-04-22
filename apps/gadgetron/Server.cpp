
#include <boost/asio.hpp>

#include <memory>

#include <Context.h>

#include "log.h"

#include "Server.h"
#include "Connection.h"
#include "connection/SocketStreamBuf.h"
#include "storage_server.h"
using namespace boost::filesystem;
using namespace Gadgetron::Server;


Server::
    Server(
        const Settings& settings
) : settings(settings) {}

[[noreturn]] void Server::serve() {


    GINFO_STREAM("Gadgetron home directory: " << settings.paths.gadgetron_home);
    GINFO_STREAM("Gadgetron working directory: " << settings.paths.working_folder);


#if(BOOST_VERSION >= 107000)
    boost::asio::io_context executor;
#else
    boost::asio::io_service executor;
#endif
    boost::asio::ip::tcp::endpoint local(boost::asio::ip::tcp::v6(), settings.port);
    boost::asio::ip::tcp::acceptor acceptor(executor, local);

    acceptor.set_option(boost::asio::socket_base::reuse_address(true));

    while(true) {
        auto socket = std::make_unique<boost::asio::ip::tcp::socket>(executor);
        acceptor.accept(*socket);

        GINFO_STREAM("Accepted connection from: " << socket->remote_endpoint().address());

        Connection::handle(settings, Gadgetron::Connection::stream_from_socket(std::move(socket)));
    }
}
