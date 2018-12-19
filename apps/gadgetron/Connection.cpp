
#include <memory>

#include "connection/ProtoConnection.h"
#include "connection/ConfigConnection.h"
#include "connection/StreamConnection.h"
#include "Connection.h"

#include "log.h"

using namespace boost::asio;
using namespace Gadgetron::Server::Connection;

namespace {

    void send_errors(std::iostream &stream) {}

    void send_close(std::iostream &stream) {
        uint16_t close = 4;
        stream.write(reinterpret_cast<char *>(&close), sizeof(close));
    }

    void handle_connection(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {

        auto config = ProtoConnection::process(*stream, paths);

        if (config) {
            auto context = ConfigConnection::process(*stream, paths);

            StreamConnection::process(*stream, context, config.get());
        }

        GDEBUG_STREAM("Sending errors.");
        send_errors(*stream);

        GDEBUG_STREAM("Sending close.");
        send_close(*stream);
    }
}

void Gadgetron::Server::Connection::handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {

    auto thread = std::thread(handle_connection, paths, std::move(stream));
    thread.detach();
}
