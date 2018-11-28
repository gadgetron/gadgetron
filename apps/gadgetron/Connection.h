#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/filesystem/path.hpp>
#include <boost/asio.hpp>

#include "Context.h"

#include "log.h"

namespace Gadgetron::Server {

    class Connection : public std::enable_shared_from_this<Connection> {

        using tcp = boost::asio::ip::tcp;

    public:
        static std::shared_ptr<Connection> create(Gadgetron::Core::Context::Paths &paths, tcp::socket &socket);

    private:
        Connection(Gadgetron::Core::Context::Paths &paths, tcp::socket &socket);

        void start();
        void process_input();
        void process_output();

        tcp::iostream stream_;

        Gadgetron::Core::Context::Paths paths_;


    };
}


#endif //GADGETRON_CONNECTION_H
