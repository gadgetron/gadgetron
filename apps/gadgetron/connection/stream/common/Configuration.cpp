#include "Configuration.h"

#include <vector>
#include <ismrmrd/ismrmrd.h>

#include "io/primitives.h"
#include "MessageID.h"

using namespace Gadgetron::Core;

namespace Gadgetron::Server::Connection::Stream {

    using Serializable = Core::variant<Config::External,Config>;

    static void send_header(std::iostream &stream, const Context::Header &header) {
        std::stringstream strstream;
        ISMRMRD::serialize(header, strstream);

        IO::write(stream, HEADER);
        IO::write_string_to_stream<uint32_t>(stream, strstream.str());
    }

    static void send_config(std::iostream &stream, const Serializable &config) {
        Core::visit([&stream](auto &config) {
                                IO::write(stream, CONFIG);
                                IO::write_string_to_stream<uint32_t>(stream, serialize_config(config));
                            },
                            config);
    }

    void Configuration::send(std::iostream &stream) const {
        send_config(stream, config);
        send_header(stream, context.header);
    }

    Configuration::Configuration(
            Core::StreamContext context,
            Config config
    ) : context(std::move(context)), config{config} {}

    Configuration::Configuration(
            Core::StreamContext context,
            Config::External config
    ) : context(std::move(context)), config{config} {}

    Configuration::Configuration(
            Core::StreamContext context,
            Config::Distributed config
    ) : Configuration(
            std::move(context),
            Config{
                config.readers,
                config.writers,
                config.stream
            }
    ) {}

    Configuration::Configuration(
            Core::StreamContext context,
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