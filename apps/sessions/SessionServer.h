#pragma once
#include <string>
#include <thread>
namespace Gadgetron::Sessions {

    std::thread start_session_server( unsigned int port );
}