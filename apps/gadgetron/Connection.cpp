
#include <memory>
#include <iostream>

#include "Context.h"

#include "connection/Core.h"
#if !(_WIN32 )
#include <unistd.h>
#endif

using namespace Gadgetron::Server::Connection;

namespace Gadgetron::Server::Connection {

#if _WIN32
    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
        auto thread = std::thread(handle_connection, std::move(stream), paths);
        thread.detach();
    }
#else
    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
        auto pid = fork();
        if (pid == 0) handle_connection(std::move(stream),paths);
    }
#endif
}
