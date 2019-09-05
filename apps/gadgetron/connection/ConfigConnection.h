#pragma once

#include "Core.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::ConfigConnection {
    void process(
        std::iostream &stream,
        const Core::Context::Paths &paths,
        const Core::Context::Args &args,
        ErrorHandler &error_handler
    );
}
