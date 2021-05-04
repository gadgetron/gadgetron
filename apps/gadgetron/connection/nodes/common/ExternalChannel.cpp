
#include "ExternalChannel.h"

#include "Serialization.h"
#include "Configuration.h"
#include "External.h"

using namespace Gadgetron::Core;

namespace Gadgetron::Server::Connection::Nodes {

    class ExternalChannel::Outbound::Closed : public ExternalChannel::Outbound {
        void push(Core::Message message) override { throw Core::ChannelClosed(); }
        void close() override {}
    };

    class ExternalChannel::Inbound::Closed : public ExternalChannel::Inbound {
        Core::Message pop() override { throw Core::ChannelClosed(); }
        void close() override {}
    };

    class ExternalChannel::Outbound::Open : public ExternalChannel::Outbound {
    public:
        explicit Open(ExternalChannel *channel) : channel(channel) {}

        void push(Core::Message message) override {
            std::lock_guard<std::mutex> guard{channel->mutex};
            channel->serialization->write(*channel->stream, std::move(message));
        }

        void close() override {
            std::lock_guard<std::mutex> guard{channel->mutex};
            channel->serialization->close(*channel->stream);
            channel->outbound = std::make_unique<Outbound::Closed>();
        }
    private:
        ExternalChannel * const channel;
    };

    class ExternalChannel::Inbound::Open : public ExternalChannel::Inbound {
    public:
        explicit Open(ExternalChannel *channel) : channel(channel) {}

        Core::Message pop() override {
            return channel->serialization->read(
                    *channel->stream,
                    [&]() {
                        auto &c_remote_errors = channel->remote_errors;
                        channel->inbound->close();

                        if (c_remote_errors.empty()) throw ChannelClosed(); else throw RemoteError(c_remote_errors);
                    },
                    [&](auto message) {
                        channel->outbound->close();
                        channel->remote_errors.push_back(message);
                    }
            );
        }

        void close() override {
            std::lock_guard<std::mutex> guard{channel->mutex};
            channel->inbound = std::make_unique<Inbound::Closed>();
        }
    private:
        ExternalChannel * const channel;
    };
}

namespace Gadgetron::Server::Connection::Nodes {

    ExternalChannel::ExternalChannel(
            std::unique_ptr<std::iostream> stream,
            std::shared_ptr<Serialization> serialization
    ) : remote_errors(),
        stream(std::move(stream)),
        serialization(std::move(serialization)) {

        inbound = std::make_unique<Inbound::Open>(this);
        outbound = std::make_unique<Outbound::Open>(this);
    }

    ExternalChannel::ExternalChannel(
            std::unique_ptr<std::iostream> stream,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration
    ) : ExternalChannel(std::move(stream), std::move(serialization)) {
        configuration->send(*this->stream);
    }

    Core::Message ExternalChannel::pop() {
        return inbound->pop();
    }

    void ExternalChannel::push_message(Core::Message message) {
        outbound->push(std::move(message));
    }

    void ExternalChannel::close() {
        outbound->close();
    }
}