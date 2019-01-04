
#include "HeaderConnection.h"

#include <map>
#include <iostream>

#include "StreamConnection.h"
#include "VoidConnection.h"
#include "Handlers.h"
#include "Config.h"

#include "io/primitives.h"
#include "Context.h"

#define CONFIG_ERROR "Received second config file. Only one allowed."

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::IO;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Handlers;

    using Header = Gadgetron::Core::Context::Header;

    class HeaderHandler : public Handler {
    public:
        explicit HeaderHandler(
                std::function<void(Header)> header_callback
        ) : header_callback(std::move(header_callback)) {}

        void handle(std::istream &stream) override {
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
        std::shared_ptr<MessageChannel> channel;
        boost::optional<Header> header;
        const Context::Paths paths;
    };

    std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(
            std::function<void()> close,
            HeaderContext &context
    ) {
        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        auto close_callback = [=, &context]() {
            context.header = boost::none;
            close();
        };

        auto header_callback = [=, &context](Header header) {
            context.header = header;
            close();
        };

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[HEADER]   = std::make_unique<HeaderHandler>(header_callback);
        handlers[QUERY]    = std::make_unique<QueryHandler>(*context.channel);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close_callback);

        return handlers;
    }
}

namespace Gadgetron::Server::Connection::HeaderConnection {

    void process(
            std::iostream &stream,
            const Core::Context::Paths &paths,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        HeaderContext context{
                std::make_shared<MessageChannel>(),
                boost::none,
                paths
        };

        std::thread input_thread = start_input_thread(
                stream,
                context.channel,
                [&](auto close) { return prepare_handlers(close, context); },
                error_handler
        );

        std::thread output_thread = start_output_thread(
                stream,
                context.channel,
                default_writers,
                error_handler
        );

        input_thread.join();
        output_thread.join();

        if (context.header) {
            StreamConnection::process(stream, Context{context.header.get(), paths}, config, error_handler);
        }
        else {
            VoidConnection::process(stream, paths, config, error_handler);
        }
    }
}