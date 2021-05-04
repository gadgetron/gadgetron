#pragma once

#include "connection/stream/common/External.h"

#include "Types.h"

namespace Gadgetron::Server::Connection::Nodes {
    std::vector<Address> discover_remote_peers();
    std::vector<Address> discover_peers();
}


