#pragma once

#include "Config.h"
#include "Builders.h"

#include "Connection.h"

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
                Context context,
                Config config
        );

    private:
        StreamConnection(std::iostream &stream, Context context, Config config);

        std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(bool &closed) override;
        std::vector<std::unique_ptr<Writer>> prepare_writers() override;

        const Config config;
        const Context context;

        Builder builder;

        std::unique_ptr<Core::Node> node;
    };
}


