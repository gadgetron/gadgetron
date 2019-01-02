#pragma once

#include <memory>

#include "Channel.h"

namespace Gadgetron::Server::Connection {
    class NodeHandler {

    public:
        virtual void process(std::shared_ptr<Core::Channel> in, std::shared_ptr<Core::Channel> out) = 0;
        virtual ~NodeHandler() = default;
    };
}


