
#include "ConfigConnection.h"

#include <map>
#include <iostream>

#include "gadgetron_config.h"

#include "Handlers.h"
#include "HeaderConnection.h"
#include "config/Config.h"

#include "io/primitives.h"
#include "Context.h"
#include "MessageID.h"
#include "Types.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::IO;
using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Handlers;

#ifdef USE_GTBABYLON
#include <GTBabylon.h>

    static std::unique_ptr<std::istream> open_and_verify_config(const std::string& filename)
    {
        auto filestream = std::ifstream(filename);
        auto config_string = std::string(std::istreambuf_iterator<char>(filestream),{});
        auto decoded =  GTBabylon::decode_message(config_string);
        return std::make_unique<std::stringstream>(decoded);
    }
#else
    static std::unique_ptr<std::istream> open_and_verify_config(const std::string& filename)
    {
        return std::make_unique<std::ifstream>(filename);
    }
#endif

namespace {

    using Header = Gadgetron::Core::StreamContext::Header;

    std::string read_filename_from_stream(std::istream &stream) {
        auto buffer = read<std::array<char,1024>>(stream);
        return std::string(buffer.data());
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
                const StreamContext::Paths &paths
        ) : ConfigHandler(callback), paths(paths) {}

        void handle(std::istream &stream, Gadgetron::Core::OutputChannel&) override {
            auto recon_name = read_filename_from_stream(stream);

            // Look up if there is an environment variable with that name
            if (getenv(recon_name.c_str())) {
                recon_name = std::string(getenv(recon_name.c_str()));
            }

            boost::filesystem::path filename = paths.gadgetron_home / GADGETRON_CONFIG_PATH / recon_name;

            GDEBUG_STREAM("Reading config file: " << filename);

            auto config_stream = open_and_verify_config(filename.string());
            handle_callback(*config_stream);
        }

    private:
        const StreamContext::Paths &paths;
    };

    class ConfigStringHandler : public ConfigHandler {
    public:
        explicit ConfigStringHandler(std::function<void(Config)> &&callback)
        : ConfigHandler(callback) {}

        void handle(std::istream &stream, Gadgetron::Core::OutputChannel& ) override {
            std::stringstream config_stream(read_string_from_stream<uint32_t>(stream));
            handle_callback(config_stream);
        }
    };

    class ConfigStreamContext {
    public:
        Gadgetron::Core::optional<Config> config;
        const StreamContext::Paths paths;
    };

    std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(
            std::function<void()> close,
            ConfigStreamContext &context
    ) {
        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        auto config_callback = [=, &context](Config config) {
            context.config = config;
            close();
        };

        handlers[FILENAME] = std::make_unique<ConfigReferenceHandler>(config_callback, context.paths);
        handlers[CONFIG]   = std::make_unique<ConfigStringHandler>(config_callback);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>("Received ISMRMRD header before config file.");
        handlers[QUERY]    = std::make_unique<QueryHandler>();
        handlers[CLOSE]    = std::make_unique<CloseHandler>(close);

        return handlers;
    }
};


namespace Gadgetron::Server::Connection::ConfigConnection {

    void process(
            std::iostream &stream,
            const Core::StreamContext::Paths &paths,
            const Core::StreamContext::Args &args,
            const Core::StreamContext::StorageAddress& sessions_address,
            ErrorHandler &error_handler
    ) {
        GINFO_STREAM("Connection state: [CONFIG]");

        ConfigStreamContext context{
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

        if (context.config) {
            HeaderConnection::process(stream, paths, args, sessions_address, context.config.value(), error_handler);
        }
    }
}