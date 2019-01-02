#pragma once

#include "Config.h"

#include "Connection.h"
#include "connection/Loader.h"

#include "Node.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {

    class StreamConnection : public Connection {
    public:
        using MessageChannel = Gadgetron::Core::MessageChannel;
        using Context = Gadgetron::Core::Context;

        static void process(
                std::iostream &stream,
                const Context &context,
                const Config &config,
                ErrorHandler &error_handler
        );

    private:
        StreamConnection(std::iostream &stream, Loader &loader);

        std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(std::function<void()> close) override;
        std::vector<std::unique_ptr<Writer>> prepare_writers() override;

        Loader &loader;
        std::unique_ptr<NodeHandler> node;
    };
}


