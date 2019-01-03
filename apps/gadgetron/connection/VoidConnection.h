#pragma once

#include "Core.h"
#include "Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::VoidConnection {
    void process(
            std::iostream &stream,
            const Core::Context::Paths &paths,
            const Config &config,
            ErrorHandler &error_handler
    );
}
