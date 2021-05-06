#pragma once

#include <memory>
#include <iostream>

#include "Context.h"
#include "Address.h"

namespace Gadgetron::Server::Connection {
    void handle(
            const Gadgetron::Core::StreamContext::Paths &paths,
            const Gadgetron::Core::StreamContext::Args &args,
            const Gadgetron::Storage::Address& storage,
            std::unique_ptr<std::iostream> stream
    );
}
