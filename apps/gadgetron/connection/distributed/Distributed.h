//
// Created by dchansen on 1/16/19.
//

#pragma once

#include <connection/NodeHandler.h>

namespace Gadgetron::Server::Distributed {
    class Distributed : public Connection::NodeHandler {

    public:
        void process(std::shared_ptr<Core::Channel> in, std::shared_ptr<Core::Channel> out) override;

    };
}



