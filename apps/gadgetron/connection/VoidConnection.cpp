#include "VoidConnection.h"

#include "Handlers.h"
#include "Writers.h"
#include "Loader.h"

#include "io/primitives.h"
#include "Reader.h"
#include "Channel.h"
#include "Context.h"
#include "RESTStorageClient.h"


namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::IO;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Writers;
    using namespace Gadgetron::Server::Connection::Handlers;

    std::vector<std::unique_ptr<Writer>> prepare_writers(std::vector<std::unique_ptr<Writer>> &writers) {
        auto ws = default_writers();

        for (auto &writer : writers) { ws.emplace_back(std::move(writer)); }

        return std::move(ws);
    }
}

namespace Gadgetron::Server::Connection::VoidConnection {

    void process(
            std::iostream &stream,
            const Core::Context::Paths &paths,
            StorageSpaces storageSpaces,
            const Config &config,
            ErrorHandler &error_handler
    ) {
        GINFO_STREAM("Connection state: [VOID]");

        // Please note the empty header initialization crime and empty storage. TODO: Fight crime.
        Core::StreamContext context{Context::Header{}, paths,storageSpaces,StreamContext::Args{}};
        Loader loader{context};


        auto ochannel = make_channel<MessageChannel>();
        auto ichannel = make_channel<MessageChannel>();

        auto node = loader.load(config.stream);
        auto writers = loader.load_writers(config);

        std::thread output_thread = start_output_thread(
                stream,
                std::move(ochannel.input),
                [&writers]() { return prepare_writers(writers); },
                error_handler
        );

        { // This.... this is not nice.
            OutputChannel removed = std::move(ichannel.output);
        }

        node->process(std::move(ichannel.input) , std::move(ochannel.output) ,error_handler);
        output_thread.join();
    }
}
