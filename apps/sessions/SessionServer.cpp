//
// Created by dch on 4/9/21.
//

#include "SessionServer.h"
#include <boost/beast.hpp>
#include <thread>
#include "RESTServer.h"

#include <nlohmann/json.hpp>
#include <boost/asio/deadline_timer.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/asio/placeholders.hpp>
#include "DB.h"

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

        template<class Read, class Send >
        void handle(Read &read, Send &&send, http::request_parser<http::empty_body>&& req) {
            json version = { {"version","0.01alpha"}};
            send(json_response(version,req.get()));
        }
    };

    struct BlobStorageEndPoint {

        std::shared_ptr<Gadgetron::Sessions::DB::DB> database;
        bool accept( const http::request<http::empty_body>& req){
            if (req.method() == http::verb::patch && req.target().starts_with(endpoint)) return true;
            return false;
        }
        template<class Read, class Send>
        void handle(Read &read, Send &&send, http::request_parser<http::empty_body>&& req) {
            beast::error_code ec;

            auto get_key = [this](const auto& req) -> std::string { auto target = req.target; target.remove_prefix(endpoint.size()); return target;};

            auto key = get_key(req.get());

            auto meta_blob = database->pending_writes[key];

            if (!meta_blob) {send(error_reponse(req.)); return;}

            http::request_parser<http::file_body> req2_parser(std::move(req));
            req2_parser.get().body().open("/tmp/testfile.bin",beast::file_mode::write,ec);
            std::cout << "Error " << ec.message() << std::endl;

            auto ec2 = read(req2_parser);
            if (ec2) std::cout <<"Error " << ec2.message() << std::endl;
            json response = {{"blob_id","testfile"}};
            send(json_response(response,req2_parser.get()));

        }

        const std::string endpoint = "/v1/blob/";

    };


}

void cleanup(const boost::system::error_code& ec, asio::deadline_timer& timer){
    std::cout << "Cleaning up all da thing" << std::endl;

    timer.expires_at(timer.expires_at()+boost::posix_time::minutes(5));
    timer.async_wait([&](const auto& ec){ cleanup(ec,timer);});
}

std::thread Gadgetron::Sessions::start_session_server( unsigned int port) {

    std::thread server_thread([=] {

        asio::io_context ioc{};
        REST::RESTEndpoints endpoints(InfoEndPoint{}, BlobStorageEndPoint{});

        asio::spawn(ioc, [&ioc,&endpoints,  port](auto yield) {
            REST::navi(ioc, tcp::endpoint(asio::ip::tcp::v4(),port), endpoints, yield);
        });

        asio::deadline_timer timer(ioc,boost::posix_time::seconds(1));
        //timer.async_wait([&](const auto& ec){ cleanup(ec,timer);});
        ioc.run();
    });
    return server_thread;
}