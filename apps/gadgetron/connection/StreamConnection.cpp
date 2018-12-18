
#include <memory>
#include <future>

#include "StreamConnection.h"

#include "Builders.h"
#include "Handlers.h"

#include "readers/Primitives.h"
#include "Reader.h"
#include "Channel.h"
#include "Context.h"


namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Handlers;

    class ReaderHandler : public Handler {
    public:
        ReaderHandler(std::unique_ptr<Reader> &&reader, std::shared_ptr<MessageChannel> channel)
                : reader(std::move(reader)), channel(std::move(channel)) {}

        void handle(std::istream &stream) override {
            channel->push_message(reader->read(stream));
        }

        std::unique_ptr<Reader> reader;
        std::shared_ptr<MessageChannel> channel;
    };

}

namespace Gadgetron::Server::Connection {

    StreamConnection::StreamConnection(
            Config config,
            Context context,
            std::unique_ptr<std::iostream> stream
    ) : config(std::move(config)), context(std::move(context)), stream(std::move(stream)) {}

    void StreamConnection::start() {

        GDEBUG_STREAM("HELLO, I'M SUPER ESCALATED!");
    }

    std::shared_ptr<StreamConnection>
    StreamConnection::create(Config config, Gadgetron::Core::Context context, std::unique_ptr<std::iostream> stream) {
        return std::make_shared<StreamConnection>(std::move(config), std::move(context), std::move(stream));
    }
}

