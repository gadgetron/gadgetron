#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/asio.hpp>

#include "Context.h"

namespace Gadgetron::Server::Connection {

    static void start(Gadgetron::Core::Context::Paths &paths, std::unique_ptr<boost::asio::ip::tcp::iostream> &stream);

}


#endif //GADGETRON_CONNECTION_H
