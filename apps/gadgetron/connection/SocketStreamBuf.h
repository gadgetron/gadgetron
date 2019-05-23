
#pragma once
#include <boost/asio/ip/tcp.hpp>
#include <iostream>

namespace Gadgetron::Connection {

    std::unique_ptr<std::iostream> stream_from_socket(std::unique_ptr<boost::asio::ip::tcp::socket> socket);
    std::unique_ptr<std::iostream> remote_stream(const std::string & host, const std::string& service);
}