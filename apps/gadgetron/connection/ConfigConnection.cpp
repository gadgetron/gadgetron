
#include "ConfigConnection.h"

#include <map>
#include <iostream>

#include "gadgetron_config.h"

#include "HeaderConnection.h"
#include "Handlers.h"
#include "Config.h"

#include "io/primitives.h"
#include "Context.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::IO;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Handlers;

    using Header = Gadgetron::Core::Context::Header;

    std::string read_filename_from_stream(std::istream &stream) {
        char buffer[1024];
        read(stream, buffer);
        return std::string(buffer);
    }

    class ConfigHandler : public Handler {
    public:
        explicit ConfigHandler(std::function<void(Config)> callback)
        : callback(std::move(callback)) {}

        void handle_callback(std::istream &config_stream) {
            callback(parse_config(config_stream));
        }

    private:
        std::function<void(Config)> callback;
    };

    class ConfigReferenceHandler : public ConfigHandler {
    public:
        ConfigReferenceHandler(
                std::function<void(Config)> &&callback,
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
        explicit ConfigStringHandler(std::function<void(Config)> &&callback)
        : ConfigHandler(callback) {}

        void handle(std::istream &stream) override {
            std::stringstream config_stream(read_string_from_stream<uint32_t>(stream));
            handle_callback(config_stream);
        }
    };

    class ConfigContext {
    public:
        std::shared_ptr<MessageChannel> channel;
        boost::optional<Config> config;
        const Context::Paths paths;
    };

    std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(
            std::function<void()> close,
            ConfigContext &context
    ) {
        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        auto config_callback = [=, &context](Config config) {
            context.config = config;
            close();
        };

        handlers[FILENAME] = std::make_unique<ConfigReferenceHandler>(config_callback, context.paths);
        handlers[CONFIG]   = std::make_unique<ConfigStringHandler>(config_callback);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>("Received ISMRMRD header before config file.");
        handlers[QUERY]    = std::make_unique<QueryHandler>(*context.channel);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close);

        return handlers;
    }
};


namespace Gadgetron::Server::Connection::ConfigConnection {

    void process(std::iostream &stream, const Core::Context::Paths &paths, ErrorHandler &error_handler) {

        ConfigContext context{
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

        if (context.config) {
            HeaderConnection::process(stream, paths, context.config.get(), error_handler);
        }
    }
}