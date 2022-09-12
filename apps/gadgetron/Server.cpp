
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
        const boost::program_options::variables_map &args,
        std::string address
) : args(args), storage_address(std::move(address)) {}

[[noreturn]] void Server::serve() {

    Gadgetron::Core::Context::Paths paths{args["home"].as<path>(), args["dir"].as<path>()};
    GINFO_STREAM("Gadgetron home directory: " << paths.gadgetron_home);
    GINFO_STREAM("Gadgetron working directory: " << paths.working_folder);

    boost::asio::io_context executor;
    boost::asio::ip::tcp::endpoint local(Info::tcp_protocol(), args["port"].as<unsigned short>());
    boost::asio::ip::tcp::acceptor acceptor(executor, local);

    acceptor.set_option(boost::asio::socket_base::reuse_address(true));

    while(true) {
        auto socket = std::make_unique<boost::asio::ip::tcp::socket>(executor);
        acceptor.accept(*socket);

        GINFO_STREAM("Accepted connection from: " << socket->remote_endpoint().address());

        Connection::handle(paths, args, storage_address, Gadgetron::Connection::stream_from_socket(std::move(socket)));
    }
}
