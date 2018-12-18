#pragma once

#include <Context.h>
#include "ProtoConnection.h"

namespace Gadgetron::Server::Connection {

    class StreamConnection : public std::enable_shared_from_this<StreamConnection> {
    private:
        StreamConnection(Core::Context context, Config config, std::unique_ptr<std::iostream> stream);

        void start();

        void process_input();


        std::unique_ptr<std::iostream> stream;

    };
}


