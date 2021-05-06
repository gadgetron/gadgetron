#ifndef GADGETRON_SERVER_H
#define GADGETRON_SERVER_H

#include <boost/asio.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options/variables_map.hpp>
#include "RESTStorageClient.h"

namespace Gadgetron::Server {

    class Server {
    public:
        Server(const boost::program_options::variables_map &args, Gadgetron::Storage::Address storage_address );

        [[noreturn]] void serve();

    private:
        const boost::program_options::variables_map &args;
        const Gadgetron::Storage::Address storage_address;
    };
}

#endif //GADGETRON_SERVER_H
