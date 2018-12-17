
#include <memory>
#include <boost/asio.hpp>

#include "connection/ProtoConnection.h"
#include "Connection.h"

using namespace boost::asio;
using namespace Gadgetron::Server::Connection;


std::shared_ptr<Connection>
Connection::create(Gadgetron::Core::Context::Paths &paths, std::unique_ptr<tcp::iostream> &stream) {

    auto connection = std::make_shared<ProtoConnection>(paths, stream);
    connection->start();

    return connection;
}
