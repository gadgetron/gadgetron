#pragma once

#include "Core.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::ConfigConnection {
    void process(
        std::iostream &stream,
        const Core::StreamContext::Paths &paths,
        const Core::StreamContext::Args &args,
        ErrorHandler &error_handler
    );
}
