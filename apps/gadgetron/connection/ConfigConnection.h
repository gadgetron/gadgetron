#pragma once

#include "Connection.h"
#include "Config.h"
#include "Channel.h"
#include "Context.h"


namespace Gadgetron::Server::Connection {

    class ConfigConnection : public Connection {
    public:
        using Context = Gadgetron::Core::Context;
        using Header = Context::Header;

        static Context process(
                std::iostream &stream,
                const Core::Context::Paths &paths
        );

    protected:
        ConfigConnection(std::iostream &stream, Gadgetron::Core::Context::Paths paths);

        std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(bool &closed) override;

        std::promise<Header> promise;

        const Gadgetron::Core::Context::Paths paths;
    };
}
