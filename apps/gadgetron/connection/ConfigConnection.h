#pragma once

#include "Connection.h"
#include "Config.h"
#include "Channel.h"
#include "Context.h"


namespace Gadgetron::Server::Connection {

    class ConfigConnection {
    public:
        using MessageChannel = Gadgetron::Core::MessageChannel;
        using Header = Gadgetron::Core::Context::Header;
        using Context = Gadgetron::Core::Context;

        static Context process(
                std::iostream &stream,
                const Core::Context::Paths &paths
        );

    private:
        ConfigConnection(Context::Paths paths, std::iostream &stream);
        ~ConfigConnection();

        void process_input();
        void process_output();

        std::promise<Header> promise;

        std::shared_ptr<MessageChannel> channel;
        std::iostream &stream;

        struct {
            std::thread input, output;
        } threads;

        const Gadgetron::Core::Context::Paths paths;
    };


}
