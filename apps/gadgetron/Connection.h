#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/asio.hpp>

#include "connection/Handlers.h"

#include "Context.h"

namespace Gadgetron::Server::Connection {
    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream);
}


#endif //GADGETRON_CONNECTION_H
