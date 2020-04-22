#pragma once

#include "Core.h"
#include "Config.h"

#include "Context.h"
#include "Connection.h"

namespace Gadgetron::Server::Connection::StreamConnection {

    void process(
            std::iostream &stream,
            const StreamContext &context,
            const Config &config,
            ErrorHandler &error_handler
    );



}


