
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

#include <nlohmann/json.hpp>

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::IO;
using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Handlers;
using json = nlohmann::json;

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

            std::ifstream file(filename.string(), std::ios::in | std::ios::binary);
            if (!file.is_open())
            {
                GDEBUG_STREAM("--> Failed to open file at path: " + filename.string());
                GDEBUG_STREAM("--> Let's check the text input ...  ");

                MessageID id = MessageID::ERROR;
                stream.read(reinterpret_cast<char*>(&id), sizeof(MessageID));

                if (id != MessageID::TEXT)
                    throw std::runtime_error("The 2nd attempt to config the chain failed ... ");

                std::string str = IO::read_string_from_stream<uint32_t>(stream);
                GDEBUG_STREAM("2nd attempt, get the config string as : " << str);

                json j = json::parse(str);
                std::string Run_This_If_Set = j["parameters"]["Run_This_If_Set"];
                std::string Select_One_To_Run = j["parameters"]["Select_One_To_Run"];

                std::string config_xml_name_from_para = Run_This_If_Set;
                if (Run_This_If_Set.empty())
                    config_xml_name_from_para = Select_One_To_Run;

                config_xml_name_from_para += ".xml";
                filename = paths.gadgetron_home / GADGETRON_CONFIG_PATH / config_xml_name_from_para;

                GDEBUG_STREAM("2nd attempt, Run_This_If_Set is " << Run_This_If_Set << " - Select_One_To_Run is " << Select_One_To_Run << " -- config file name is " << filename.string());
            }

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