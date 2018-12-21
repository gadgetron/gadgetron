#pragma once

#include "Config.h"
#include "Builders.h"


#include "Node.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {

    class StreamConnection  {
    public:
        using MessageChannel = Gadgetron::Core::MessageChannel;
        using Context = Gadgetron::Core::Context;

    StreamConnection(std::iostream &stream, Context context, Config config);
    ~StreamConnection();
       void process();
    private:

        std::thread output_thread;
        std::thread stream_thread;
        std::iostream& stream;
        std::shared_ptr<MessageChannel> input_channel = std::make_shared<MessageChannel>();
        std::shared_ptr<MessageChannel> output_channel = std::make_shared<MessageChannel>();

        const Config config;
        const Context context;

        Builder builder;

    };
}


