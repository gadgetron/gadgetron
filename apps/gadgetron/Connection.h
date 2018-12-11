#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/filesystem/path.hpp>
#include <boost/asio.hpp>

#include "Context.h"

#include "log.h"

namespace Gadgetron::Server {

    class Connection {
        using tcp = boost::asio::ip::tcp;

    public:
        static std::shared_ptr<Connection> create(Gadgetron::Core::Context::Paths &paths, std::unique_ptr<tcp::iostream> &stream);
        virtual ~Connection() = default;
    };
}


#endif //GADGETRON_CONNECTION_H
