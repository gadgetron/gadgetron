#pragma once

#include "Config.h"
#include "Builders.h"

#include "Node.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {

    class StreamConnection : public std::enable_shared_from_this<StreamConnection> {
    public:
        using MessageChannel = Gadgetron::Core::MessageChannel;
        using Context = Gadgetron::Core::Context;

        static void process(
                std::iostream &stream,
                Context context,
                Config config
        );

    private:

        StreamConnection(Gadgetron::Core::Context context, Config config, std::iostream &stream);
        ~StreamConnection();

        void process_input();
        void process_output();

        const Config config;
        const Context context;

        struct {
            std::shared_ptr<MessageChannel> input, output;
        } channels;

        struct {
            std::thread input, output;
        } threads;

        std::iostream &stream;

        Builder builder;

        std::unique_ptr<Core::Node> node;
    };
}


