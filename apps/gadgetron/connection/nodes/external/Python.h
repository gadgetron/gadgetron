#pragma once

#ifdef __clang__
   #define  unary_function  __unary_function
#endif

#include <boost/process.hpp>

#include "connection/config/Config.h"

#include "Context.h"

namespace Gadgetron::Server::Connection::Nodes {
    boost::process::child start_python_module(
        const Config::Execute &,
        unsigned short port,
        const Gadgetron::Core::StreamContext &
    );
    bool python_available() noexcept;
}

