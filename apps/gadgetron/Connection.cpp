
#include <memory>
#include <iostream>

#include "Context.h"

#include "connection/Core.h"

using namespace Gadgetron::Server::Connection;

namespace Gadgetron::Server::Connection {

    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
        auto thread = std::thread(handle_connection, std::move(stream), paths);
        thread.detach();
    }
}
