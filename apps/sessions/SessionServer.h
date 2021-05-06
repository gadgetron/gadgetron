#pragma once
#include <string>
#include <thread>
#include <boost/filesystem.hpp>
#include <boost/asio/io_context.hpp>

namespace Gadgetron::Storage {

    class SessionServer {
    public:
        SessionServer(unsigned short port, const boost::filesystem::path& database_folder, const boost::filesystem::path& blob_folder);
        ~SessionServer();
        unsigned short port();

    private:
        std::thread server_thread;
        boost::asio::io_context ioContext;
        unsigned short bound_port;

    };

}