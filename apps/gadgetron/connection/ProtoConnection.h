#pragma once

#include "Connection.h"
#include "Channel.h"
#include "Context.h"
#include "Builders.h"

namespace Gadgetron::Server::Connection {

    class ProtoConnection : public Connection, public std::enable_shared_from_this<ProtoConnection> {
    public:
        using tcp = boost::asio::ip::tcp;
        using MessageChannel = Gadgetron::Core::MessageChannel;

        ProtoConnection(Gadgetron::Core::Context::Paths paths, std::unique_ptr<tcp::iostream> &stream)
                : stream(std::move(stream))
                , paths(std::move(paths)) {}

        void start();
        void process_input();
        void process_output();

        void escalate(Config config);

        const Gadgetron::Core::Context::Paths paths;
        std::shared_ptr<MessageChannel> channel;
        std::unique_ptr<tcp::iostream> stream;
    };

}

