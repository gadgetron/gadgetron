//
// Created by dchansen on 9/10/19.
//
#include <random>

#ifdef __clang__
    #define unary_function  __unary_function
#endif

#include <boost/asio.hpp>
#include <gtest/gtest.h>

#include "../connection/SocketStreamBuf.h"
#include "log.h"


namespace ba = boost::asio;
using tcp    = boost::asio::ip::tcp;
using namespace Gadgetron;


class SocketTest : public ::testing::Test {

public:
    SocketTest() : ::testing::Test() {
        auto endpoint = tcp::endpoint(tcp::v6(), 0);
        acceptor = std::make_unique<tcp::acceptor>(tcp::acceptor(ios, endpoint));

        auto port    = acceptor->local_endpoint().port();
        auto socketF = std::async([&]() {
            tcp::socket socket{ ios };
            acceptor->accept(socket);
            return socket;
        });

        socketstream = Connection::remote_stream("localhost", std::to_string(port));

        server_socket = std::make_unique<tcp::socket>(socketF.get());
    }


    ba::io_service ios{};
    std::unique_ptr<tcp::acceptor> acceptor;
    std::unique_ptr<std::iostream> socketstream;
    std::unique_ptr<tcp::socket> server_socket;

};

TEST_F(SocketTest, read_test) {

    auto data = std::vector<char>(19,42);
    ba::write(*server_socket,ba::buffer(data.data(),data.size()));

    auto data2 = std::vector<char>(19);

    socketstream->read(data2.data(),data2.size());

    ASSERT_EQ(data,data2);
}

TEST_F(SocketTest, write_test) {

    auto data = std::vector<char>(1u << 22,42);

    auto thread = std::thread([&](){socketstream->write(data.data(),data.size());});

    auto data2 = std::vector<char>(data.size());
    ba::read(*server_socket,ba::buffer(data2.data(),data2.size()));

    ASSERT_EQ(data,data2);
    thread.join();
}


TEST_F(SocketTest, stringstream_test) {
    const std::string name = "Albatros";
    std::stringstream sstream;
    sstream << name;
    *socketstream << sstream.rdbuf();



    std::vector<char> ref(name.size());
    ba::read(*server_socket,ba::buffer(ref.data(),ref.size()));

    std::string name2(ref.begin(),ref.end());
    ASSERT_EQ(name,name2);
}



TEST_F(SocketTest, nulltest) {
    auto data =  std::vector<char>(1u << 22,0);
    std::mt19937_64 engine;
    std::uniform_int_distribution<int> distribution(0);
    for (auto& d : data) d = distribution(engine);

    GINFO_STREAM(std::to_string(data.size()) << std::endl);
    std::stringstream sstream;
    sstream.write(data.data(),data.size());

    auto thread = std::thread([&](){   *socketstream << sstream.rdbuf(); *socketstream << sstream.rdbuf();});



    std::vector<char> ref(data.size());
    ba::read(*server_socket,ba::buffer(ref.data(),ref.size()));

    ASSERT_EQ(ref,data);
    thread.join();
}
