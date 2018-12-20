
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
#include "Connection_common.h"


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

    void StreamConnection::process() {

        auto factory = [&](auto& closed ) {

            std::function<void()> close_callback = [&]() {
                closed = true;
                input_channel->close();
            };

            std::map<uint16_t, std::unique_ptr<Handler>> handlers{};
            handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
            handlers[CONFIG] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
            handlers[HEADER] = std::make_unique<ErrorProducingHandler>(HEADER_ERROR);

            handlers[QUERY] = std::make_unique<QueryHandler>(input_channel);
            handlers[CLOSE] = std::make_unique<CloseHandler>(close_callback);

            for (auto &reader : builder.build_readers(config.readers)) {
                handlers[reader.first] = std::make_unique<ReaderHandler>(std::move(reader.second), input_channel);
            }

            return handlers;
        };
    }


    StreamConnection::StreamConnection(std::iostream &stream, Gadgetron::Core::Context context, Config config)
    : stream(stream), context(context), config(config), builder(context.paths) {

        output_thread = std::thread(
                [&](auto writers){ handle_output(output_channel,stream,std::move(writers));},
                builder.build_writers(this->config.writers)
                );

        stream_thread = std::thread(
                [this](){
                    auto node = builder.build_stream(this->config.stream, this->context);
                    node->process(this->input_channel,this->output_channel);
                });
    }

    StreamConnection::~StreamConnection() {
        input_channel->close();
        stream_thread.join();
        output_thread.join();

    }
}

