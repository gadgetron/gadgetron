#pragma once

#include "../../distributed/RemoteChannel.h"

namespace Gadgetron::Server::Connection::Stream {

    struct Local   { };
    struct Address { std::string ip, port; };

    using  Worker = Core::variant<Address, Local>;

    std::vector<Worker> get_remote_workers();
    std::vector<Worker> get_workers();
}


