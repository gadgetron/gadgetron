
#include <memory>

#include "StreamConnection.h"

#include "Handlers.h"
#include "Writers.h"
#include "Loader.h"

#include "readers/Primitives.h"
#include "Reader.h"
#include "Channel.h"
#include "Context.h"


namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Writers;
    using namespace Gadgetron::Server::Connection::Handlers;

    class ReaderHandler : public Handler {
    public:
        ReaderHandler(std::unique_ptr<Reader> &&reader, std::shared_ptr<MessageChannel> channel)
                : reader(std::move(reader)), channel(std::move(channel)) {}

        void handle(std::istream &stream) override {
            channel->push_message(reader->read(stream));
        }

        std::unique_ptr<Reader> reader;
        std::shared_ptr<MessageChannel> channel;
    };

}

#define CONFIG_ERROR "Received second config file. Only one config allowed."
#define HEADER_ERROR "Received second ISMRMRD header. Only one allowed."

namespace Gadgetron::Server::Connection {

    std::map<uint16_t, std::unique_ptr<Connection::Handler>> StreamConnection::prepare_handlers(std::function<void()> close) {

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>(HEADER_ERROR);

        handlers[QUERY]    = std::make_unique<QueryHandler>(*channels.input);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close);

        for (auto &reader : loader.readers()) {
            handlers[reader.first] = std::make_unique<ReaderHandler>(std::move(reader.second), channels.input);
        }

        return handlers;
    }

    std::vector<std::unique_ptr<Writer>> StreamConnection::prepare_writers() {
        return loader.writers();
    }

    StreamConnection::StreamConnection(std::iostream &stream, Loader &loader)
            : Connection(stream), loader(loader) {

        channels.input = std::make_shared<MessageChannel>();
        channels.output = std::make_shared<MessageChannel>();

        node = loader.stream();
    }

    void StreamConnection::process(
            std::iostream &stream,
            const Context &context,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        Loader loader{error_handler, context, config};
        StreamConnection connection{stream, loader};

        std::thread input_thread = error_handler.run(
                "Connection Input Thread",
                [&]() { connection.process_input(); }
        );

        std::thread output_thread = error_handler.run(
                "Connection Output Thread",
                [&]() { connection.process_output(); }
        );

        connection.node->process(connection.channels.input, connection.channels.output);

        input_thread.join();
        output_thread.join();
    }
}

