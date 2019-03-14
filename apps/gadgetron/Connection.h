#pragma once

#include <memory>
#include <iostream>

#include "Context.h"

namespace Gadgetron::Server::Connection {
    void handle(
            const Gadgetron::Core::Context::Paths &paths,
            const Gadgetron::Core::Context::Args &args,
            std::unique_ptr<std::iostream> stream
    );
}
