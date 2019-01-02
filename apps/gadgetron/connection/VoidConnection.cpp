#include "VoidConnection.h"

#include "Handlers.h"
#include "Writers.h"
#include "Loader.h"

#include "readers/Primitives.h"
#include "Reader.h"
#include "Channel.h"
#include "Context.h"


namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Writers;
    using namespace Gadgetron::Server::Connection::Handlers;

}

namespace Gadgetron::Server::Connection {

    std::map<uint16_t, std::unique_ptr<Connection::Handler>> VoidConnection::prepare_handlers(std::function<void()> close) {
        throw std::runtime_error("VoidConnection has no handlers - it doesn't have an input loop.");
    }


    std::vector<std::unique_ptr<Writer>> VoidConnection::prepare_writers() {
        return loader.writers();
    }

    VoidConnection::VoidConnection(std::iostream &stream, Loader &loader)
          : Connection(stream), loader(loader) {

        channels.input = std::make_shared<MessageChannel>();
        channels.output = std::make_shared<MessageChannel>();

        node = loader.stream();
    }

    void VoidConnection::process(
            std::iostream &stream,
            const Context::Paths &paths,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        ISMRMRD::IsmrmrdHeader header{};
        Context context{header, paths};

        Loader loader{error_handler, context, config};
        VoidConnection connection{stream, loader};

        std::thread output_thread = error_handler.run(
                "Connection Output Thread",
                [&]() { connection.process_output(); }
        );

        connection.node->process(connection.channels.input, connection.channels.output);
        output_thread.join();
    }
}
