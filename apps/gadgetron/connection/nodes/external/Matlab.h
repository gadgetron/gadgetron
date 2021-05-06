#pragma once

#include <boost/process.hpp>

#include "connection/config/Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::Nodes {
    boost::process::child start_matlab_module(const Config::Execute &, unsigned short port, const Gadgetron::Core::Context &);
    bool matlab_available() noexcept;
}