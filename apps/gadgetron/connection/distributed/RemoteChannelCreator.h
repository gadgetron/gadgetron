//
// Created by dchansen on 1/16/19.
//

#pragma once

#include "distributed/ChannelCreator.h"
namespace Gadgetron::Server::Distributed {

    class RemoteChannelCreator : public Core::Distributed::ChannelCreator {
    public:
        Core::OutputChannel &create() override;


    private:

    };

}


