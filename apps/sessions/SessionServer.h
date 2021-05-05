#pragma once
#include <string>
#include <thread>
#include <filesystem>
#include <boost/asio/io_context.hpp>

namespace Gadgetron::Sessions {

    class SessionServer {
    public:
        SessionServer(unsigned short port, const std::filesystem::path& database_folder, const std::filesystem::path& blob_folder);
        ~SessionServer();
        unsigned short port();

    private:
        std::thread server_thread;
        boost::asio::io_context ioContext;
        unsigned short bound_port;

    };

}