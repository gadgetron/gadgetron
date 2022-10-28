#include "Connection.h"
#include <iostream>
#include <memory>

#include "Context.h"

#include "connection/Core.h"
#if !(_WIN32)
#include <cstdlib>
#include <unistd.h>
#include <sys/wait.h>
#endif

using namespace Gadgetron::Server::Connection;

namespace Gadgetron::Server::Connection {

#if _WIN32 || !NDEBUG || GADGETRON_DISABLE_FORK || __clang__

    void handle(
            const Gadgetron::Core::StreamContext::Paths& paths,
            const Gadgetron::Core::StreamContext::Args& args,
            const std::string& storage_address,
            std::unique_ptr<std::iostream> stream
    ) {
        auto thread = std::thread(handle_connection, std::move(stream), paths, args, storage_address);
        thread.detach();
    }

#else

    void handle(
            const Gadgetron::Core::StreamContext::Paths& paths,
            const Gadgetron::Core::StreamContext::Args& args,
            const Gadgetron::Core::StreamContext::StorageAddress& storage_address,
            std::unique_ptr<std::iostream> stream
    ) {
        auto pid = fork();
        if (pid == 0) {
            handle_connection(std::move(stream), paths, args, storage_address);
            std::quick_exit(0);
        }
        auto listen_for_close = [](auto pid) {int status; waitpid(pid,&status,0);};
        std::thread t(listen_for_close,pid);
        t.detach();
    }

#endif
}
