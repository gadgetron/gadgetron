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

    class VoidContext {
    public:
        struct {
            std::shared_ptr<MessageChannel> input, output;
        } channels;
        Loader loader;
    };

    std::vector<std::unique_ptr<Writer>> prepare_writers(Loader &loader) {
        auto writers = loader.writers();
        auto defaults = default_writers();

        for (auto &writer : defaults) { writers.emplace_back(std::move(writer)); }

        return std::move(writers);
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
        VoidContext ctx{
                {std::make_shared<MessageChannel>(), std::make_shared<MessageChannel>()},
                Loader{error_handler, Context{Context::Header{}, paths}, config}
        };

        auto node = ctx.loader.stream();

        std::thread output_thread = start_output_thread(
                stream,
                ctx.channels.output,
                [&]() { return prepare_writers(ctx.loader); },
                error_handler
        );

        ctx.channels.input->close();
        node->process(ctx.channels.input, ctx.channels.output);
        output_thread.join();
    }
}
