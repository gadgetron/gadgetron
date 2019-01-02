#pragma once

#include "Connection.h"
#include "Config.h"
#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection {

    class ConfigConnection : public Connection {
    public:
        static void process(
                std::iostream &stream,
                const Core::Context::Paths &paths,
                ErrorHandler &error_handler
        );

    protected:
        ConfigConnection(std::iostream &stream, Gadgetron::Core::Context::Paths paths);

        std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(std::function<void()>) override;

        std::promise<boost::optional<Config>> promise;
        const Gadgetron::Core::Context::Paths paths;
    };

}

