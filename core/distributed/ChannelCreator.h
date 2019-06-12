#pragma once

#include <memory>
#include "../Channel.h"

namespace Gadgetron::Core::Distributed {
    class ChannelCreator {
    public:
        virtual OutputChannel create() = 0;
        virtual ~ChannelCreator() = default;
    };
}



