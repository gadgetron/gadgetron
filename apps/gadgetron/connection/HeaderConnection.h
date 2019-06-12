#pragma once

#include "Core.h"
#include "Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::HeaderConnection {
    void process(
            std::iostream &stream,
            const Core::Context::Paths &paths,
            const Core::Context::Args &args,
            const Config &config,
            ErrorHandler &error_handler
    );
}
