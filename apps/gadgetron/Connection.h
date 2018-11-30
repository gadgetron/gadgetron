#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/filesystem/path.hpp>
#include <boost/asio.hpp>

#include "Stream.h"
#include "Context.h"
#include "Channel.h"

#include "Config.h"

#include "log.h"

namespace Gadgetron::Server {

    class Connection {
        using tcp = boost::asio::ip::tcp;

    public:
        static std::shared_ptr<Connection> create(Gadgetron::Core::Context::Paths &paths, tcp::socket &socket);
        virtual ~Connection() = default;
    };
}


#endif //GADGETRON_CONNECTION_H
