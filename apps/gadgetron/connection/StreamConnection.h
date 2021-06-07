#pragma once

#include "Core.h"
#include "config/Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::StreamConnection {
    void process(
            std::iostream &stream,
            const Core::StreamContext &context,
            const Config &config,
            ErrorHandler &error_handler
    );
}


