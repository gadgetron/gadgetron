#pragma once

#include "Core.h"
#include "config/Config.h"

#include "Context.h"

using namespace Gadgetron::Core;

namespace Gadgetron::Server::Connection::StreamConnection {
    void process(
            std::iostream &stream,
            const StreamContext &context,
            const Config &config,
            ErrorHandler &error_handler
    );
}


