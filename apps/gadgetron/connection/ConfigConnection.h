#pragma once

#include "Core.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::ConfigConnection {
    void process(
        std::iostream &stream,
        const Settings& settings,
        ErrorHandler &error_handler
    );
}
