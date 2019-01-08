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

        for (auto &writer : writers) { writers.emplace_back(std::move(writer)); }

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
        // Please note the header crime. TODO: Fight crime.
        Context context{Context::Header{}, paths};
        Loader loader{error_handler, context, config};

        struct {
            std::shared_ptr<MessageChannel> input = std::make_shared<MessageChannel>();
            std::shared_ptr<MessageChannel> output = std::make_shared<MessageChannel>();
        } channels;

        auto node = loader.stream();
        auto writers = loader.writers();

        std::thread output_thread = start_output_thread(
                stream,
                channels.output,
                [&]() { return prepare_writers(writers); },
                error_handler
        );

        channels.input->close();
        node->process(channels.input, channels.output);
        output_thread.join();
    }
}
