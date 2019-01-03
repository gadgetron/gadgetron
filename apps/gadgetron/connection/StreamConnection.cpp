
#include <memory>

#include "StreamConnection.h"

#include "Handlers.h"
#include "Writers.h"
#include "Loader.h"

#include "readers/Primitives.h"
#include "Reader.h"
#include "Channel.h"
#include "Context.h"

#define CONFIG_ERROR "Received second config file. Only one allowed."
#define HEADER_ERROR "Received second ISMRMRD header. Only one allowed."

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

    class StreamContext {
    public:
        struct {
            std::shared_ptr<MessageChannel> input, output;
        } channels;

        Loader loader;
    };

    std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(
            std::function<void()> close,
            StreamContext &context
    ) {
        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>(HEADER_ERROR);
        handlers[QUERY]    = std::make_unique<QueryHandler>(*context.channels.input);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close);

        for (auto &reader : context.loader.readers()) {
            handlers[reader.first] = std::make_unique<ReaderHandler>(
                    std::move(reader.second),
                    context.channels.input
            );
        }

        return handlers;
    }

    std::vector<std::unique_ptr<Writer>> prepare_writers(Loader &loader) {
        auto writers = loader.writers();
        auto defaults = default_writers();

        for (auto &writer : defaults) { writers.emplace_back(std::move(writer)); }

        return std::move(writers);
    }
}


namespace Gadgetron::Server::Connection::StreamConnection {

    void process(
            std::iostream &stream,
            const Core::Context &context,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        StreamContext ctx{
                {std::make_shared<MessageChannel>(), std::make_shared<MessageChannel>()},
                Loader{error_handler, context, config}
        };

        auto node = ctx.loader.stream();

        std::thread input_thread = start_input_thread(
                stream,
                ctx.channels.input,
                [&](auto close) { return prepare_handlers(close, ctx); },
                error_handler
        );

        std::thread output_thread = start_output_thread(
                stream,
                ctx.channels.output,
                [&]() { return prepare_writers(ctx.loader); },
                error_handler
        );

        node->process(ctx.channels.input, ctx.channels.output);

        input_thread.join();
        output_thread.join();
    }
}
