#pragma once

#include <boost/process.hpp>

#include "connection/config/Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::Nodes {
    boost::process::child start_julia_module(
        const Config::Execute &,
        unsigned short port,
        const Gadgetron::Core::StreamContext &
    );
    bool julia_available() noexcept;
}

