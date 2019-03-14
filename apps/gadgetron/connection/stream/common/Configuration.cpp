#include "Configuration.h"

#include <vector>
#include <ismrmrd/ismrmrd.h>

#include "io/primitives.h"
#include "MessageID.h"

using namespace Gadgetron::Core;

namespace Gadgetron::Server::Connection::Stream {

    void send_header(std::iostream &stream, const Context::Header &header) {
        std::stringstream strstream;
        ISMRMRD::serialize(header, strstream);

        IO::write(stream, HEADER);
        IO::write_string_to_stream<uint32_t>(stream, strstream.str());
    }

    void send_config(std::iostream &stream, const Configuration::Serializable &serializable) {
        serializable.write(stream);
    }

    void Configuration::send(std::iostream &stream) const {
        send_config(stream, *config);
        send_header(stream, context.header);
    }

    class ConfigSerializable : public Configuration::Serializable {
    public:
        const Config config;

        explicit ConfigSerializable(Config config) : config(std::move(config)) {}

        void write(std::iostream &stream) const override {
            IO::write(stream, CONFIG);
            IO::write_string_to_stream<uint32_t>(stream, serialize_config(config));
        };
    };

    Configuration::Configuration(
            Core::Context context,
            Config config
    ) : context(std::move(context)), config(std::make_unique<ConfigSerializable>(config)) {}

    Configuration::Configuration(
            Core::Context context,
            Config::PureDistributed config
    ) : Configuration(
            std::move(context),
            Config{
                config.readers,
                config.writers,
                Config::Stream {
                    "PureStream",
                    std::vector<Config::Node>(config.stream.gadgets.begin(), config.stream.gadgets.end())
                }
            }
        ) {}
}