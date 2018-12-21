#pragma once

#include "Connection.h"
#include "Config.h"
#include "Channel.h"
#include "Context.h"

#include "Connection_common.h"
namespace Gadgetron::Server::Connection {

    using ConfigConnection = BasicConnection<Core::Context::Header>;

    template<> Core::Context::Header ConfigConnection::process();

}
