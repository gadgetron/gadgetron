//
// Created by dch on 4/9/21.
//

#include "SessionServer.h"
#include <cstdint>
#include <boost/beast.hpp>
using namespace Gadgetron::Sessions;
namespace beast = boost::beast;         // from <boost/beast.hpp>
namespace http = beast::http;           // from <boost/beast/http.hpp>
namespace net = boost::asio;            // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>


void test() {
    http::verb::post
};