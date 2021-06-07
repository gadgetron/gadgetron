
#include "HeaderConnection.h"

#include <map>
#include <iostream>

#include "Handlers.h"
#include "StreamConnection.h"
#include "VoidConnection.h"
#include "config/Config.h"

#include "io/primitives.h"
#include "Context.h"
#include "MessageID.h"

#define CONFIG_ERROR "Received second config file. Only one allowed."

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::IO;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Handlers;

    using Header = Gadgetron::Core::StreamContext::Header;

    class HeaderHandler : public Handler {
    public:
        explicit HeaderHandler(
                std::function<void(Header)> header_callback
        ) : header_callback(std::move(header_callback)) {}

        void handle(std::istream &stream, OutputChannel&) override {
            std::string raw_header(read_string_from_stream<uint32_t>(stream));

            ISMRMRD::IsmrmrdHeader header{};
            ISMRMRD::deserialize(raw_header.c_str(), header);

            header_callback(header);
        }

    private:
        std::function<void(Header)> header_callback;
    };

    class HeaderContext {
    public:
        Gadgetron::Core::optional<Header> header;
        const StreamContext::Paths paths;
    };

    std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(
            std::function<void()> close,
            HeaderContext &context
    ) {
        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        auto header_callback = [=, &context](Header header) {
            context.header = header;
            close();
        };

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[HEADER]   = std::make_unique<HeaderHandler>(header_callback);
        handlers[QUERY]    = std::make_unique<QueryHandler>();
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close);

        return handlers;
    }
}

namespace Gadgetron::Server::Connection::HeaderConnection {

    void process(
            std::iostream &stream,
            const Core::StreamContext::Paths &paths,
            const Core::StreamContext::Args &args,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        GINFO_STREAM("Connection state: [HEADER]");

        HeaderContext context{
                Core::none,
                paths
        };

        auto channel = make_channel<MessageChannel>();

        std::thread input_thread = start_input_thread(
                stream,
                std::move(channel.output),
                [&](auto close) { return prepare_handlers(close, context); },
                error_handler
        );

        std::thread output_thread = start_output_thread(
                stream,
                std::move(channel.input),
                default_writers,
                error_handler
        );

        input_thread.join();
        output_thread.join();

        if (context.header) {
            StreamConnection::process(stream, StreamContext{context.header.value(), paths, args}, config, error_handler);
        }
        else {
            VoidConnection::process(stream, paths, config, error_handler);
        }
    }
}