#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/asio.hpp>

#include "Context.h"

namespace Gadgetron::Server::Connection {

    class Connection {
        using tcp = boost::asio::ip::tcp;

    public:
        virtual ~Connection() = default;
        static std::shared_ptr<Connection> create(Gadgetron::Core::Context::Paths &paths, std::unique_ptr<tcp::iostream> &stream);
    };
}


#endif //GADGETRON_CONNECTION_H
