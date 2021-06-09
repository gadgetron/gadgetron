#pragma once

#include <boost/asio.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options/variables_map.hpp>

namespace Gadgetron::Server {

    class Server {
    public:
        Server(const boost::program_options::variables_map &args, std::string storage_address);

        [[noreturn]] void serve();

    private:
        const boost::program_options::variables_map &args;
        const std::string storage_address;
    };
}
