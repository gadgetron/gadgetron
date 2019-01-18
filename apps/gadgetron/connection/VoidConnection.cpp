#include "VoidConnection.h"

#include "Handlers.h"
#include "Writers.h"
#include "Loader.h"

#include "io/primitives.h"
#include "Reader.h"
#include "Channel.h"
#include "Context.h"


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
            const Config &config,
            ErrorHandler &error_handler
    ) {
        GINFO_STREAM("Connection state: [VOID]");

        // Please note the empty header initialization crime. TODO: Fight crime.
        Context context{Context::Header{}, paths};
        Loader loader{context};

        struct {
            std::shared_ptr<MessageChannel> input = std::make_shared<MessageChannel>();
            std::shared_ptr<MessageChannel> output = std::make_shared<MessageChannel>();
        } channels;

        auto node = loader.load(config.stream);
        auto writers = loader.load_writers(config);

        std::thread output_thread = start_output_thread(
                stream,
                channels.output,
                [&]() { return prepare_writers(writers); },
                error_handler
        );

        channels.input->close();
        node->process(channels.input, channels.output, error_handler);
        output_thread.join();
    }
}
