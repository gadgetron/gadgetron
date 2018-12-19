
#include "ProtoConnection.h"

#include <typeindex>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/dll/shared_library.hpp>
#include <boost/dll.hpp>
#include <boost/range/algorithm/transform.hpp>

#include "log.h"
#include "gadgetron_config.h"

#include "readers/Primitives.h"
#include "Response.h"
#include "Reader.h"
#include "Writer.h"

#include "Server.h"
#include "Config.h"

#include "Writers.h"
#include "Handlers.h"
#include "StreamConnection.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Writers;
    using namespace Gadgetron::Server::Connection::Handlers;

    using Header = Gadgetron::Core::Context::Header;

    std::string read_filename_from_stream(std::istream &stream) {
        char buffer[1024];
        read_into(stream, buffer);
        return std::string(buffer);
    }

    class ConfigHandler : public Handler {
    public:
        explicit ConfigHandler(std::function<void(Config)> &callback)
        : callback(callback) {}

        void handle_callback(std::istream &config_stream) {
            callback(parse_config(config_stream));
        }

    private:
        std::function<void(Config)> &callback;
    };

    class ConfigReferenceHandler : public ConfigHandler {
    public:
        ConfigReferenceHandler(
                std::function<void(Config)> &callback,
                const Context::Paths &paths
        ) : ConfigHandler(callback), paths(paths) {}

        void handle(std::istream &stream) override {
            boost::filesystem::path filename = paths.gadgetron_home / GADGETRON_CONFIG_PATH / read_filename_from_stream(stream);

            GDEBUG_STREAM("Reading config file: " << filename << std::endl);

            std::ifstream config_stream(filename.string());
            handle_callback(config_stream);
        }

    private:
        const Context::Paths &paths;
    };

    class ConfigStringHandler : public ConfigHandler {
    public:
        explicit ConfigStringHandler(std::function<void(Config)> &callback)
        : ConfigHandler(callback) {}

        void handle(std::istream &stream) override {
            std::stringstream config_stream(read_string_from_stream<uint32_t>(stream));
            handle_callback(config_stream);
        }
    };


};

namespace Gadgetron::Server::Connection {

    void ProtoConnection::process_input() {

        bool closed = false;

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        std::function<void()> close_callback = [&]() {

            closed = true;
            promise.set_value(boost::none);
        };

        std::function<void(Config)> config_callback = [&](Config config) {

            closed = true;
            promise.set_value(boost::make_optional(config));
        };

        handlers[FILENAME] = std::make_unique<ConfigReferenceHandler>(config_callback, paths);
        handlers[CONFIG]   = std::make_unique<ConfigStringHandler>(config_callback);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>("Received ISMRMRD header before config file.");
        handlers[QUERY]    = std::make_unique<QueryHandler>(channel);
        handlers[CLOSE]    = std::make_unique<CallbackHandler>(close_callback);

        while (!closed) {
            auto id = read_t<uint16_t>(stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(stream);
        }
    }

    void ProtoConnection::process_output() {

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

    ProtoConnection::ProtoConnection(Context::Paths paths, std::iostream &stream)
    : stream(stream), paths(std::move(paths)), channel(std::make_shared<MessageChannel>()) {}

    boost::optional<Config> ProtoConnection::process(std::iostream &stream, const Context::Paths &paths) {

        ProtoConnection connection{paths, stream};

        connection.threads.input  = std::thread([&]() {
            connection.process_input();
        });

        connection.threads.output = std::thread([&]() {
            connection.process_output();
        });

        auto future = connection.promise.get_future();
        return future.get();
    }

    ProtoConnection::~ProtoConnection() {

        // Terminate the input thread somehow? Sabotage the stream?
        channel->close();

        threads.input.join();
        threads.output.join();
    }
}