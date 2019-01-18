
#include <memory>

#include "StreamConnection.h"

#include "Handlers.h"
#include "Writers.h"
#include "Loader.h"

#include "io/primitives.h"
#include "Reader.h"
#include "Channel.h"
#include "Context.h"
#include "MessageID.h"

#define CONFIG_ERROR "Received second config file. Only one allowed."
#define HEADER_ERROR "Received second ISMRMRD header. Only one allowed."

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::IO;
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

    std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(
            std::function<void()> close,
            std::shared_ptr<MessageChannel> input,
            std::map<uint16_t, std::unique_ptr<Reader>> &readers
    ) {
        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>(HEADER_ERROR);
        handlers[QUERY]    = std::make_unique<QueryHandler>(*input);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close);

        for (auto &pair : readers) {
            handlers[pair.first] = std::make_unique<ReaderHandler>(std::move(pair.second), input);
        }

        return handlers;
    }

    std::vector<std::unique_ptr<Writer>> prepare_writers(std::vector<std::unique_ptr<Writer>> &writers) {
        auto ws = default_writers();

        for (auto &writer : writers) { ws.emplace_back(std::move(writer)); }

        return ws;
    }
}


namespace Gadgetron::Server::Connection::StreamConnection {

    void process(
            std::iostream &stream,
            const Core::Context &context,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        GINFO_STREAM("Connection state: [STREAM]");

        Loader loader{context};

        struct {
            std::shared_ptr<MessageChannel> input = std::make_shared<MessageChannel>();
            std::shared_ptr<MessageChannel> output = std::make_shared<MessageChannel>();
        } channels;

        auto node = loader.load(config.stream);
        auto readers = loader.load_readers(config);
        auto writers = loader.load_writers(config);

        std::thread input_thread = start_input_thread(
                stream,
                channels.input,
                [&](auto close) { return prepare_handlers(close, channels.input, readers); },
                error_handler
        );

        std::thread output_thread = start_output_thread(
                stream,
                channels.output,
                [&]() { return prepare_writers(writers); },
                error_handler
        );

        node->process(channels.input, channels.output, error_handler);

        input_thread.join();
        output_thread.join();
    }
}
