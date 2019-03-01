
#include "ExternalChannel.h"

#include "io/primitives.h"
#include "MessageID.h"

using namespace Gadgetron::Core;

namespace {
    std::string make_error_message(const std::vector<std::string>& errors) {
        std::stringstream error_maker;
        error_maker << "Errors received: " << std::endl;
        for (auto& error : errors) {
            error_maker << error << std::endl;
        }
        return error_maker.str();
    }
}

namespace Gadgetron::Server::Connection::Stream {

    RemoteError::RemoteError(std::vector<std::string> messages) : std::runtime_error(make_error_message(messages)), messages(messages) {}

    class ExternalChannel::Outbound {
    public:
        virtual ~Outbound() = default;
        virtual void push(Core::Message message) = 0;
        virtual void close() = 0;

        class Open; class Closed;
    };

    class ExternalChannel::Inbound {
    public:
        virtual ~Inbound() = default;
        virtual Core::Message pop() = 0;
        virtual void close() = 0;

        class Open; class Closed;
    };

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
            channel->serialization.write(*channel->stream, std::move(message));
        }

        void close() override {
            std::lock_guard guard{channel->mutex};
            channel->outbound = std::make_unique<Outbound::Closed>();
            IO::write(*channel->stream, CLOSE);
        }
    private:
        ExternalChannel *channel;
    };

    class ExternalChannel::Inbound::Open : public ExternalChannel::Inbound {
    public:
        explicit Open(ExternalChannel *channel) : channel(channel) {}

        Core::Message pop() override {
            return channel->serialization.read(
                    *channel->stream,
                    [&]() {
                        channel->inbound->close();
                        if (channel->remote_errors.empty()) throw ChannelClosed();
                        throw RemoteError(channel->remote_errors);
                    },
                    [&](auto message) {
                        channel->outbound->close();
                        channel->remote_errors.push_back(message);
                    }
            );
        }

        void close() override {
            std::lock_guard guard{channel->mutex};
            channel->inbound = std::make_unique<Inbound::Closed>();
        }
    private:
        ExternalChannel *channel;
    };
}

namespace Gadgetron::Server::Connection::Stream {

    Serialization::Serialization(
            Serialization::Readers readers,
            Serialization::Writers writers
    ) : readers(std::move(readers)), writers(std::move(writers)) {}

    void Serialization::write(std::iostream &stream, Core::Message message) const {

        auto writer = std::find_if(
                writers.begin(),
                writers.end(),
                [&](auto& writer) { return writer->accepts(message); }
        );

        if (writer == writers.end())
            throw std::runtime_error("Could not find appropriate writer for message in external channel.");

        (*writer)->write(stream, std::move(message));
    }

    Core::Message Serialization::read(
            std::iostream &stream,
            std::function<void()> on_close,
            std::function<void(std::string message)> on_error
    ) const {

        auto id = IO::read<uint16_t>(stream);
        auto illegal_message = [&](auto &) {
            throw std::runtime_error("Received illegal message id from external peer: " + std::to_string(id));
        };

        std::map<uint16_t, std::function<void(std::iostream &)>> handlers{
                {FILENAME,  illegal_message},
                {CONFIG,    illegal_message},
                {HEADER,    illegal_message},
                {CLOSE,     [&](auto &) { on_close(); }},
                {TEXT,      illegal_message},
                {QUERY,     illegal_message},
                {RESPONSE,  illegal_message},
                {ERROR,     [&](auto &stream) { on_error(IO::read_string_from_stream<uint64_t>(stream)); }}
        };

        for (; handlers.count(id); id = IO::read<uint16_t>(stream)) handlers.at(id)(stream);

        return readers.at(id)->read(stream);
    }

    Configuration::Configuration(
            Core::Context context,
            Config::External external
    ) : context(std::move(context)), external(std::move(external)) {}

    ExternalChannel::ExternalChannel(
            std::unique_ptr<std::iostream> stream,
            Serialization serialization,
            Configuration configuration
    ) : stream(std::move(stream)),
        serialization(std::move(serialization)),
        configuration(std::move(configuration)) {

        inbound = std::make_unique<Inbound::Open>(this);
        outbound = std::make_unique<Outbound::Open>(this);
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