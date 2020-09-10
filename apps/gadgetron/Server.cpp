
#include <boost/asio.hpp>


#include <Context.h>

#include "log.h"

#include "Server.h"
#include "Connection.h"
#include "connection/SocketStreamBuf.h"
#include "system_info.h"

using namespace boost::filesystem;
using namespace Gadgetron::Server;


Server::Server(
        const boost::program_options::variables_map &args
) : args(args) {}

void Server::serve() {

    Gadgetron::Core::Context::Paths paths{args["home"].as<path>(), args["dir"].as<path>()};
    GINFO_STREAM("Gadgetron home directory: " << paths.gadgetron_home);
    GINFO_STREAM("Gadgetron working directory: " << paths.working_folder);

#if(BOOST_VERSION >= 107000)
    boost::asio::io_context executor;
#else
    boost::asio::io_service executor;
#endif
    boost::asio::ip::tcp::endpoint local(Info::tcp_protocol(), args["port"].as<unsigned short>());
    boost::asio::ip::tcp::acceptor acceptor(executor, local);

    acceptor.set_option(boost::asio::socket_base::reuse_address(true));

    while(true) {
        auto socket = std::make_unique<boost::asio::ip::tcp::socket>(executor);
        acceptor.accept(*socket);

        GINFO_STREAM("Accepted connection from: " << socket->remote_endpoint().address());

        Connection::handle(paths, args, Gadgetron::Connection::stream_from_socket(std::move(socket)));
    }
}
