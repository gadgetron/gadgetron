#pragma once

#include "connection/Config.h"

#include "Context.h"
#include "Types.h"

namespace Gadgetron::Server::Connection::Stream {

    class Configuration {
    public:
        const Core::StreamContext context;

        void send(std::iostream &stream) const;

        Configuration(Core::StreamContext context, Config config);
        Configuration(Core::StreamContext context, Config::External config);
        Configuration(Core::StreamContext context, Config::Distributed config);
        Configuration(Core::StreamContext context, Config::PureDistributed config);

    private:
        Core::variant<Config::External,Config> config;
    };
}