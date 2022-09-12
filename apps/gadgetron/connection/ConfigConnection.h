#pragma once

#include "Core.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::ConfigConnection {
    void process(
        std::iostream &stream,
        const Core::StreamContext::Paths &paths,
        const Core::StreamContext::Args &args,
        const Core::StreamContext::StorageAddress& session_address,
        ErrorHandler &error_handler
    );
}
