//
// Created by dch on 4/9/21.
//

#include "SessionServer.h"
#include <boost/beast.hpp>
#include <thread>
#include "RESTServer.h"

#include <nlohmann/json.hpp>
using namespace Gadgetron::Sessions;
namespace beast = boost::beast;         // from <boost/beast.hpp>
namespace http = beast::http;           // from <boost/beast/http.hpp>
namespace net = boost::asio;            // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>
namespace asio = boost::asio;

using json = nlohmann::json;

namespace {
    template<class REQ>
    auto json_response(const json& object, const REQ& req){
        http::response<http::string_body> response(http::status::ok, req.version());
        response.set(http::field::server, BOOST_BEAST_VERSION_STRING);
        response.set(http::field::content_type, "application/json");
        response.keep_alive(req.keep_alive());
        response.body() = object.dump();
        response.prepare_payload();
        return response;
    }
    struct InfoEndPoint {

        bool accept(const http::request<http::empty_body> &req) {
            if (req.method() == http::verb::get && req.target() == "/v1/info") return true;
            return false;
        }

        template<class SyncReadStream, class Send, class Buffer>
        void handle(SyncReadStream &stream, Send &&send, Buffer &buffer, asio::yield_context yield,
                    http::request_parser<http::empty_body>&& req) {
            json version = { {"version","0.01alpha"}};
            send(json_response(version,req.get()));
        }
    };

    struct BlobStorageEndPoint {

        bool accept( const http::request<http::empty_body>& req){
            if (req.method() == http::verb::patch && req.target() == "/v1/blob") return true;
            return false;
        }
        template<class SyncReadStream, class Send, class Buffer>
        void handle(SyncReadStream &stream, Send &&send, Buffer &buffer, asio::yield_context yield,
                    http::request_parser<http::empty_body>&& req) {
            beast::error_code ec;

            http::request_parser<http::file_body> req2_parser(std::move(req));
            req2_parser.get().body().open("/tmp/testfile.bin",beast::file_mode::write,ec);
            std::cout << "Error " << ec.message() << std::endl;
            http::async_read(stream,buffer,req2_parser,yield[ec]);

            json response = {{"blob_id","testfile"}};
            send(json_response(response,req2_parser.get()));


        }

    };


}

std::thread Gadgetron::Sessions::start_session_server( unsigned int port) {

    std::thread server_thread([=] {

        asio::io_context ioc{};
        REST::RESTEndpoints endpoints(InfoEndPoint{}, BlobStorageEndPoint{});

        asio::spawn(ioc, [&ioc,&endpoints,  port](auto yield) {
            REST::navi(ioc, tcp::endpoint(asio::ip::tcp::v4(),port), endpoints, yield);
        });
        ioc.run();
    });
    return server_thread;
}