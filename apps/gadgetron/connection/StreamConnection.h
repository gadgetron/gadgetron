#pragma once

#include "Core.h"
#include "Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::StreamConnection {
    void process(
            std::iostream &stream,
            const Core::Context &context,
            const Config &config,
            ErrorHandler &error_handler
    );
}


