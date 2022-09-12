#pragma once

#include "Core.h"
#include "config/Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::HeaderConnection {
    void process(
            std::iostream &stream,
            const Core::StreamContext::Paths &paths,
            const Core::StreamContext::Args &args,
            const Core::StreamContext::StorageAddress& address,
            const Config &config,
            ErrorHandler &error_handler
    );
}
