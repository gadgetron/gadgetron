#pragma once

#include <Channel.h>
#include "Connection.h"
#include "Config.h"
#include "Channel.h"
#include "Context.h"
#include "Connection_common.h"
namespace Gadgetron::Server::Connection {

    using ProtoConnection = Connection<boost::optional<Config>>;

    template<> boost::optional<Config> ProtoConnection::process();

}

