#pragma once

#include <boost/asio.hpp>

#include "connection/Handlers.h"



#include "Writer.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {
    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream);
}


