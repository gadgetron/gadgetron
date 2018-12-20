
#include <readers/Primitives.h>
#include "ConfigConnection.h"

#include "Context.h"

#include "Handlers.h"
#include "Writers.h"
#include "Config.h"


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

#define CONFIG_ERROR "Received second config file. Only one config allowed."

namespace Gadgetron::Server::Connection {

    std::map<uint16_t, std::unique_ptr<Connection::Handler>> ConfigConnection::prepare_handlers(bool &closed) {

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        std::function<void(Header)> deliver = [&](Header header) {
            closed = true;
            channels.input->close();
            promise.set_value(header);
        };

        std::function<void()> close_callback = [&, deliver]() {

            // This crime is done to support void chains - chains with no input.
            // Examples include dependency query gadgets. TODO: Fight crime.
            Header header{};
            deliver(header);
        };

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[QUERY]    = std::make_unique<QueryHandler>(channels.output);

        handlers[HEADER]   = std::make_unique<HeaderHandler>(deliver);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close_callback);

        return handlers;
    }

    Context ConfigConnection::process(
            std::iostream &stream,
            const Core::Context::Paths &paths,
            ErrorHandler &error_handler
    ) {
        ConfigConnection connection{stream, paths};

        connection.start(error_handler);
        connection.join();

        auto future = connection.promise.get_future();
        return Context{ future.get(), paths };
    }

    ConfigConnection::ConfigConnection(std::iostream &stream, Context::Paths paths)
    : Connection(stream), paths(std::move(paths)) {
        channels.input = channels.output = std::make_shared<MessageChannel>();
    }
}