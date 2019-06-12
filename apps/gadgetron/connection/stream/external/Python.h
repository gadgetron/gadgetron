#pragma once

#include "connection/Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::Stream {
    void start_python_module(const Config::Execute &, unsigned short port, const Gadgetron::Core::Context &);
}

