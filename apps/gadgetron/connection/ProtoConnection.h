#pragma once

#include "Connection.h"
#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection {

    class ProtoConnection : public std::enable_shared_from_this<ProtoConnection> {
    public:
        static std::shared_ptr<ProtoConnection> create(
                Gadgetron::Core::Context::Paths paths,
                std::unique_ptr<std::iostream> stream
        );
        void start();

        ProtoConnection(Gadgetron::Core::Context::Paths paths, std::unique_ptr<std::iostream> stream);

    private:
        using MessageChannel = Gadgetron::Core::MessageChannel;

        void process_input();
        void process_output();

        const Gadgetron::Core::Context::Paths paths;

        std::shared_ptr<MessageChannel> channel;
        std::unique_ptr<std::iostream> stream;

        struct {
            std::thread input, output;
        } threads;
    };

}

