#pragma once

#include <memory>
#include <iostream>

#include "Context.h"

namespace Gadgetron::Server::Connection {
    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream);
}
