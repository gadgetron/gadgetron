
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
#include "Connection_common.h"

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


}

namespace Gadgetron::Server::Connection {

     template<> boost::optional<Config> ProtoConnection::process()  {


        boost::optional<Config> config = boost::none;

         auto factory = [&](auto& closed){

            std::function<void(Config)> config_callback = [&config,&closed](Config input_config) {
                config = input_config;
                closed = true;
            };
             std::map<uint16_t, std::unique_ptr<Handler>> handlers{};
             handlers[FILENAME] = std::make_unique<ConfigReferenceHandler>(config_callback, paths);
             handlers[CONFIG] = std::make_unique<ConfigStringHandler>(config_callback);
             handlers[HEADER] = std::make_unique<ErrorProducingHandler>("Received ISMRMRD header before config file.");
             handlers[QUERY] = std::make_unique<QueryHandler>(channel);

             return handlers;
         };

         handle_input(stream,factory);

        return config;

    }

    template class BasicConnection<boost::optional<Config>>;


}