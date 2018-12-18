#pragma once

#include "Config.h"

#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {

    class StreamConnection : public std::enable_shared_from_this<StreamConnection> {
    public:
        static std::shared_ptr<StreamConnection> create(
                Config config,
                Gadgetron::Core::Context context,
                std::unique_ptr<std::iostream> stream
        );
        void start();

        StreamConnection(Config config, Gadgetron::Core::Context context, std::unique_ptr<std::iostream> stream);

    private:
        using MessageChannel = Gadgetron::Core::MessageChannel;

        void process_input();
        void process_output();

        const Config config;
        const Gadgetron::Core::Context context;

        struct {
            std::shared_ptr<MessageChannel> input, output;
        } channels;

        struct {
            std::thread input, output;
        } threads;

        std::unique_ptr<std::iostream> stream;
    };
}


