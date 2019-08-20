#pragma once

#include "connection/Config.h"

#include "Context.h"
#include "Types.h"

namespace Gadgetron::Server::Connection::Stream {

    class Configuration {
    public:
        const Core::Context context;

        void send(std::iostream &stream) const;

        Configuration(Core::Context context, Config config);
        Configuration(Core::Context context, Config::External config);
        Configuration(Core::Context context, Config::Distributed config);
        Configuration(Core::Context context, Config::PureDistributed config);

    private:
        Core::variant<Config::External,Config> config;
    };
}