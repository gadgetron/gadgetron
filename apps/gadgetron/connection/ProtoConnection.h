#pragma once

#include "Connection.h"
#include "Config.h"
#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection {

    class ProtoConnection : public Connection {
    public:
        static boost::optional<Config> process(
                std::iostream &stream,
                const Core::Context::Paths &paths
        );

    protected:
        ProtoConnection(std::iostream &stream, Gadgetron::Core::Context::Paths paths);

        std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(bool &closed) override;

        std::promise<boost::optional<Config>> promise;

        const Gadgetron::Core::Context::Paths paths;
    };

}

