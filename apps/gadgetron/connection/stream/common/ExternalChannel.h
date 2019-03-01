#pragma once

#include <map>

#include "connection/Config.h"

#include "Context.h"
#include "Channel.h"
#include "Reader.h"
#include "Writer.h"

namespace Gadgetron::Server::Connection::Stream {

    class RemoteError : public std::runtime_error {
    public:
        explicit RemoteError(std::vector<std::string> messages);

    private:
        const std::vector<std::string> messages;
    };

    class Serialization {
    public:
        using Readers = std::map<uint16_t,std::unique_ptr<Core::Reader>>;
        using Writers = std::vector<std::unique_ptr<Core::Writer>>;

        const Readers readers;
        const Writers writers;

        Serialization(Readers readers, Writers writers);

        void write(std::iostream &stream, Core::Message message) const;
        Core::Message read(
                std::iostream &stream,
                std::function<void()> on_close,
                std::function<void(std::string message)> on_error
        ) const;
    };

    class Configuration {
    public:
        const Core::Context context;
        const Config::External external;

        Configuration(Core::Context context, Config::External external);
    };

    class ExternalChannel {
    public:
        ExternalChannel(
                std::unique_ptr<std::iostream> stream,
                Serialization serialization,
                Configuration configuration
        );

        Core::Message pop();

        void push_message(Core::Message message);
        void close();

    private:
        std::mutex mutex;
        std::unique_ptr<std::iostream> stream;
        std::vector<std::string> remote_errors;

        const Serialization serialization;
        const Configuration configuration;

        class Outbound; class Inbound;
        friend Outbound; friend Inbound;

        std::unique_ptr<Outbound> outbound;
        std::unique_ptr<Inbound>  inbound;
    };
}