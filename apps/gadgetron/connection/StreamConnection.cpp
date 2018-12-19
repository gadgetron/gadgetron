
#include <memory>
#include <future>

#include "StreamConnection.h"

#include "Builders.h"
#include "Handlers.h"
#include "Writers.h"

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

    void StreamConnection::process_input() {

        bool closed = false;

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        std::function<void()> close_callback = [&]() {
            closed = true;
            channels.input->close();
        };

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>(HEADER_ERROR);

        handlers[QUERY]    = std::make_unique<QueryHandler>(channels.output);
        handlers[CLOSE]    = std::make_unique<CallbackHandler>(close_callback);

        for (auto &reader : builder.build_readers(config.readers)) {
            handlers[reader.first] = std::make_unique<ReaderHandler>(std::move(reader.second), channels.input);
        }

        while (!closed) {
            auto id = read_t<uint16_t>(stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(stream);
        }
    }

    void StreamConnection::process_output() {

        std::vector<std::unique_ptr<Writer>> writers = builder.build_writers(config.writers);

        writers.push_back(std::make_unique<TextWriter>());
        writers.push_back(std::make_unique<ResponseWriter>());

        InputChannel<Message>& output = *channels.output;
        for (auto message : output) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                   [&](auto &writer) { return writer->accepts(*message); }
            );

            (*writer)->write(stream, std::move(message));
        }
    }


    void StreamConnection::process(std::iostream &stream, Context context, Config config) {

        StreamConnection connection(std::move(context), std::move(config), stream);

        connection.threads.input  = std::thread([&]() {
            connection.process_input();
        });

        connection.threads.output = std::thread([&]() {
            connection.process_output();
        });

        connection.node->process(connection.channels.input, connection.channels.output);

        // Connection destructor joins input and output threads.
    }

    StreamConnection::StreamConnection(Gadgetron::Core::Context context, Config config, std::iostream &stream)
    : context(context), config(config), stream(stream), channels {
        std::make_shared<MessageChannel>(),
        std::make_shared<MessageChannel>()
    }, builder(context.paths) {
        node = builder.build_stream(config.stream, context);
    }

    StreamConnection::~StreamConnection() {

        channels.output->close();

        threads.input.join();
        threads.output.join();
    }
}

