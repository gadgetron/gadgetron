
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
                std::function<void(Header)> &header_callback
        ) : header_callback(header_callback) {}

        void handle(std::istream &stream) override {
            std::string raw_header(read_string_from_stream<uint32_t>(stream));

            ISMRMRD::IsmrmrdHeader header{};
            ISMRMRD::deserialize(raw_header.c_str(), header);

            header_callback(header);
        }

    private:
        std::function<void(Header)> &header_callback;
    };
}

#define CONFIG_ERROR "Received second config file. Only one config allowed."

namespace Gadgetron::Server::Connection {

    void ConfigConnection::process_input() {

        GDEBUG_STREAM("Input thread running.");

        bool closed = false;

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        std::function<void()> close_callback = [&]() {

            // This crime is done to support void chains - chains with no input.
            // Examples include dependency query gadgets. TODO: Fight crime.
            Header header{};

            closed = true;
            promise.set_value(header);
        };

        std::function<void(Header)> header_callback = [&](Header header) {

            closed = true;
            promise.set_value(header);
        };

        handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
        handlers[QUERY]    = std::make_unique<QueryHandler>(channel);

        handlers[HEADER]   = std::make_unique<HeaderHandler>(header_callback);
        handlers[CLOSE]    = std::make_unique<CallbackHandler>(close_callback);

        while (!closed) {
            auto id = read_t<uint16_t>(stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(stream);
        }
    }

    void ConfigConnection::process_output() {

        std::vector<std::unique_ptr<Writer>> writers{};

        writers.push_back(std::make_unique<TextWriter>());
        writers.push_back(std::make_unique<ResponseWriter>());

        InputChannel<Message>& input = *channel;
        for (auto message : input) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                                       [&](auto &writer) { return writer->accepts(*message); }
            );

            (*writer)->write(stream, std::move(message));
        }
    }

    Context ConfigConnection::process(std::iostream &stream, const Core::Context::Paths &paths) {

        ConfigConnection connection{paths, stream};

        connection.threads.input  = std::thread([&]() {
            connection.process_input();
        });

        connection.threads.output = std::thread([&]() {
            connection.process_output();
        });

        auto future = connection.promise.get_future();
        return Context{ future.get(), paths };
    }

    ConfigConnection::ConfigConnection(
            Gadgetron::Core::Context::Paths paths,
            std::iostream &stream
    ) : paths(std::move(paths)), stream(stream), channel(std::make_shared<MessageChannel>()) {}

    ConfigConnection::~ConfigConnection() {

        // Terminate the input thread somehow? Sabotage the stream?
        channel->close();

        threads.input.join();
        threads.output.join();
    }
}