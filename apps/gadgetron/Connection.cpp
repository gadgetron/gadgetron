
#include <memory>

#include "connection/ProtoConnection.h"
#include "Connection.h"

using namespace boost::asio;
using namespace Gadgetron::Server::Connection;

void Gadgetron::Server::Connection::handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
    auto connection = ProtoConnection::create(paths, std::move(stream));
    connection->start();
}
