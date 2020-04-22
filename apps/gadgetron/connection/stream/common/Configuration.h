#pragma once

#include "connection/Config.h"
#include "connection/StreamConnection.h"

#include "Context.h"
#include "Types.h"

namespace Gadgetron::Server::Connection::Stream {

    class Configuration {
    public:
        const StreamContext context;

        void send(std::iostream &stream) const;

        Configuration(StreamContext context, Config config);
        Configuration(StreamContext context, Config::External config);
        Configuration(StreamContext context, Config::Distributed config);
        Configuration(StreamContext context, Config::PureDistributed config);

    private:
        Core::variant<Config::External,Config> config;
    };
}