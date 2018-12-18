
#include <memory>
#include <boost/asio.hpp>

#include "connection/ProtoConnection.h"
#include "Connection.h"

using namespace boost::asio;
using namespace Gadgetron::Server::Connection;

using tcp = boost::asio::ip::tcp;

void Gadgetron::Server::Connection::start(Gadgetron::Core::Context::Paths &paths, std::unique_ptr<tcp::iostream> &stream) {
    auto connection = std::make_shared<ProtoConnection>(paths, stream);
    connection->start();
}
