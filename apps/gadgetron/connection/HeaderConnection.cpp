
#include "HeaderConnection.h"

#include "readers/Primitives.h"
#include "Context.h"

#include "Handlers.h"
#include "Writers.h"
#include "Config.h"

#include "StreamConnection.h"
#include "VoidConnection.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Writers;
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
}

#define CONFIG_ERROR "Received second config file. Only one allowed."

namespace Gadgetron::Server::Connection {

    std::map<uint16_t, std::unique_ptr<Connection::Handler>> HeaderConnection::prepare_handlers(std::function<void()> close) {

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        auto deliver = [=](boost::optional<Header> header) {
            close();
            promise.set_value(header);
        };

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[QUERY]    = std::make_unique<QueryHandler>(*channels.input);

        handlers[HEADER]   = std::make_unique<HeaderHandler>([=](Header header) {
            deliver(boost::make_optional(header));
        });

        handlers[CLOSE]    = std::make_unique<CloseHandler>([=]() {
            deliver(boost::none);
        });

        return handlers;
    }

    void HeaderConnection::process(
            std::iostream &stream,
            const Context::Paths &paths,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        HeaderConnection connection{stream, paths};

        std::thread input_thread = error_handler.run(
                "Connection Input Thread",
                [&]() { connection.process_input(); }
        );

        std::thread output_thread = error_handler.run(
                "Connection Output Thread",
                [&]() { connection.process_output(); }
        );

        input_thread.join();
        output_thread.join();

        auto future = connection.promise.get_future();
        auto header = future.get();

        if (header) {
            Context context{header.get(), paths};
            StreamConnection::process(stream, context, config, error_handler);
        }
        else {
            VoidConnection::process(stream, paths, config, error_handler);
        }
    }

    HeaderConnection::HeaderConnection(std::iostream &stream, Context::Paths paths)
    : Connection(stream), paths(std::move(paths)) {
        channels.input = channels.output = std::make_shared<MessageChannel>();
    }
}