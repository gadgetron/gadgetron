#pragma once

#include "ProtoConnection.h"

namespace Gadgetron::Server::Connection {

    class StreamConnection : public std::enable_shared_from_this<StreamConnection> {
    private:
        StreamConnection(Context context, Config config) : builder(config) {

        }

        Builder builder;
    };
}


