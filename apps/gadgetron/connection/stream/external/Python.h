#pragma once

#include "Context.h"

namespace Gadgetron::Server::Connection::Stream {
    void start_python_module(std::string module, unsigned short port, const Gadgetron::Core::Context &);
}

