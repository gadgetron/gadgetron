#pragma once

#include <map>
#include <list>

#include "connection/Config.h"

#include "Context.h"
#include "Channel.h"
#include "Reader.h"
#include "Writer.h"

#include "Serialization.h"
#include "Configuration.h"

namespace Gadgetron::Server::Connection::Stream {

    class ExternalChannel {
    public:
        ExternalChannel(
                std::unique_ptr<std::iostream> stream,
                std::shared_ptr<Serialization> serialization
        );

        ExternalChannel(
                std::unique_ptr<std::iostream> stream,
                std::shared_ptr<Serialization> serialization,
                std::shared_ptr<Configuration> configuration
        );

        Core::Message pop();
        void push_message(Core::Message message);
        void close();

    private:
        std::unique_ptr<std::iostream> stream;
        std::list<std::string> remote_errors;

        std::shared_ptr<Serialization> serialization;

        class Outbound {
        public:
            virtual ~Outbound() = default;
            virtual void push(Core::Message message) = 0;
            virtual void close() = 0;

            class Open; class Closed;
        };

        class Inbound {
        public:
            virtual ~Inbound() = default;
            virtual Core::Message pop() = 0;
            virtual void close() = 0;

            class Open; class Closed;
        };

        friend Outbound; friend Inbound;

        std::mutex mutex;
        std::unique_ptr<Outbound> outbound;
        std::unique_ptr<Inbound>  inbound;
    };
}