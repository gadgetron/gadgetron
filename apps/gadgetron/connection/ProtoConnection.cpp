
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

    class HeaderHandler : public Handler {
    public:
        explicit HeaderHandler(
            std::function<void(Config, Header)> &header_callback,
            Config config
        ) : header_callback(header_callback), config(std::move(config)) {}

        void handle(std::istream &stream) override {
            std::string raw_header(read_string_from_stream<uint32_t>(stream));

            ISMRMRD::IsmrmrdHeader header{};

            if (raw_header.empty()) {
                GWARN_STREAM(
                        "Received empty ISMRMRD header from client. " <<
                        "This is deprecated, and only allowed for backwards compatibility. " <<
                        "This will not be allowed in the future."
                );
            }
            else {
                ISMRMRD::deserialize(raw_header.c_str(), header);
            }

            header_callback(config, header);
        }

    private:
        std::function<void(Config, Header)> &header_callback;
        Config config;
    };

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

    void ProtoConnection::start() {
        auto self = shared_from_this();

        threads.input = std::thread ([=]() {
            self->process_input();
        });

        threads.output = std::thread([=]() {
            self->process_output();
        });

        threads.input.detach();
    };


    void ProtoConnection::process_input() {

        GDEBUG_STREAM("Input thread running.");

        bool closed = false;

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        std::function<void(Config, Header)> header_callback = [&](Config config, Header header) {

            // I should call the close handler with a post-close callback.
            // Close handler should do the joining.

            closed = true;
            channel->close();
            threads.output.join();

            Context context{header, paths};

            auto connection = StreamConnection::create(config, context, std::move(stream));
            connection->start();
        };

        std::function<void(Config)> config_callback = [&](Config config) {
            std::string error = "Received second config file. Only one config allowed.";

            handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(error);
            handlers[CONFIG]   = std::make_unique<ErrorProducingHandler>(error);
            handlers[HEADER]   = std::make_unique<HeaderHandler>(header_callback, config);
        };

        handlers[FILENAME] = std::make_unique<ConfigReferenceHandler>(config_callback, paths);
        handlers[CONFIG]   = std::make_unique<ConfigStringHandler>(config_callback);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>("Received ISMRMRD header before config file.");
        handlers[CLOSE]    = std::make_unique<CloseHandler>(closed);
        handlers[QUERY]    = std::make_unique<QueryHandler>(channel);

        while (!closed) {
            auto id = read_t<uint16_t>(*stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(*stream);
        }

//        channel->close();
    }

    void ProtoConnection::process_output() {
        GDEBUG_STREAM("Output thread running.");

        std::vector<std::unique_ptr<Writer>> writers{};

        writers.push_back(std::make_unique<TextWriter>());
        writers.push_back(std::make_unique<ResponseWriter>());

        InputChannel<Message>& input = *channel;
        for (auto message : input) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                    [&](auto &writer) { return writer->accepts(*message); }
            );

            (*writer)->write(*stream, std::move(message));
        }
    }

    std::shared_ptr<ProtoConnection>
    ProtoConnection::create(Gadgetron::Core::Context::Paths paths, std::unique_ptr<std::iostream> stream) {
        return std::make_shared<ProtoConnection>(paths, std::move(stream));
    }

    ProtoConnection::ProtoConnection(Context::Paths paths, std::unique_ptr<std::iostream> stream)
    : stream(std::move(stream)), paths(std::move(paths)), channel(std::make_shared<MessageChannel>()) {}
}