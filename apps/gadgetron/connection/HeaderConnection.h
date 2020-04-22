#pragma once

#include "Core.h"
#include "Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::HeaderConnection {
    void process(
            std::iostream &stream,
            const Settings& settings,
            const Config &config,
            ErrorHandler &error_handler
    );
}
