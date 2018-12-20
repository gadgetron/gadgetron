
#include <memory>

#include "connection/ProtoConnection.h"
#include "connection/ConfigConnection.h"
#include "connection/StreamConnection.h"
#include "connection/Writers.h"

#include "Connection.h"


#include "readers/Primitives.h"
#include "log.h"

using namespace boost::asio;

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Readers;

using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Writers;

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

        send_errors(*stream);

        send_close(*stream);
    }
}

namespace Gadgetron::Server::Connection {

    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
        auto thread = std::thread(handle_connection, paths, std::move(stream));
        thread.detach();
    }

    Connection::Connection(std::iostream &stream) : stream(stream) {};

    void Connection::process_input() {

        bool closed = false;
        auto handlers = prepare_handlers(closed);

        while (!closed) {
            auto id = read_t<uint16_t>(stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(stream);
        }
    }

    void Connection::process_output() {

        auto writers = prepare_writers();

        writers.push_back(std::make_unique<TextWriter>());
        writers.push_back(std::make_unique<ResponseWriter>());

        InputChannel<Message> &output = *channels.output;
        for (auto message : output) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                   [&](auto &writer) { return writer->accepts(*message); }
            );

            (*writer)->write(stream, std::move(message));
        }
    }

    std::vector<std::unique_ptr<Writer>> Connection::prepare_writers() {
        return std::vector<std::unique_ptr<Writer>>();
    }

    void Connection::start() {
        threads.input  = std::thread([&]() {
            process_input();
        });

        threads.output = std::thread([&]() {
            process_output();
        });
    }

    void Connection::join() {
        threads.input.join();
        threads.output.join();
    }
}

