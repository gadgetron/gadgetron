#pragma once

#include "Core.h"
#include "Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::HeaderConnection {
    void process(
            std::iostream &stream,
            const Core::StreamContext::Paths &paths,
            const Core::StreamContext::Args &args,
            const Storage::Address& address,
            const Config &config,
            ErrorHandler &error_handler
    );
}
