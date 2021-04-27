//
// Created by dch on 4/26/21.
//

#include "RESTStorageClient.h"
#include <ismrmrd/xml.h>
#include <boost/beast/core.hpp>
#include <boost/beast/http.hpp>
#include <boost/beast/version.hpp>
#include <boost/asio/connect.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <nlohmann/json.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace beast = boost::beast;     // from <boost/beast.hpp>
namespace http = beast::http;       // from <boost/beast/http.hpp>
namespace net = boost::asio;        // from <boost/asio.hpp>
using tcp = net::ip::tcp;
using json = nlohmann::json;

namespace {
    using namespace Gadgetron::Core;

    constexpr int http_version = 11;
    auto make_json_request(http::verb method, const std::string& host, const std::string& target, const json& j){

        http::request<http::string_body> req{method,target,http_version};
        req.set(http::field::host,host);
        req.set(http::field::user_agent, BOOST_BEAST_VERSION_STRING);
        req.set(http::field::content_type,"application/json");
        req.body() = j.dump();
        req.prepare_payload();
        return req;
    }
    auto make_empty_request(http::verb method, const std::string& host, const std::string& target ){

        http::request<http::empty_body> req{method,target,http_version};
        req.set(http::field::host,host);
        req.set(http::field::user_agent, BOOST_BEAST_VERSION_STRING);
        return req;
    }
    auto  connect_stream(net::io_context& ioc, const std::string& host, const std::string& port){
        tcp::resolver resolver(ioc);
        beast::tcp_stream  stream(ioc);
        auto const results = resolver.resolve(host,port);
        stream.connect(results);
        return stream;
    }

    json get_content(const std::string& server_address, const std::string& port, const std::string& group, const std::string& subject, const std::string& key){

        net::io_context ioc;
        auto stream = connect_stream(ioc,server_address,port);
        auto json_body = json{{"storagespace", group},{"subject",subject},{"key",key}};
        auto req = make_json_request(http::verb::get,server_address,"/v1/data",json_body );
        http::write(stream,req);

        beast::flat_buffer buffer;
        http::response<http::string_body> response;
        http::read(stream,buffer, response);
        stream.socket().shutdown(tcp::socket::shutdown_both);
        return json::parse(response.body());
    }
    std::vector<char> fetch_data(const std::string& server_address, const std::string port, const std::string& path){
        net::io_context ioc;
        auto stream = connect_stream(ioc,server_address,port);
        auto req = make_empty_request(http::verb::get, server_address, path);
        http::write(stream,req);
        beast::flat_buffer buf;

        http::response_parser<http::vector_body<char>> response_parser;
        response_parser.body_limit(128ull*1024ull*1024ull*1024ull); //We support files up to 128GB. For now.
        http::read(stream,buf,response_parser);
        stream.socket().shutdown(tcp::socket::shutdown_both);
        return std::move(response_parser.get().body());
    }

    json store_request(const std::string& server_address, const std::string& port, const std::string& group, const std::string& subject, const std::string& key, const boost::posix_time::time_duration& duration){

        net::io_context ioc;
        auto stream = connect_stream(ioc,server_address,port);
        auto json_body = json{{"storagespace", group},{"subject",subject},{"key",key},{"duration", boost::posix_time::to_iso_string(duration)}};
        auto req = make_json_request(http::verb::post,server_address,"/v1/data",json_body );
        http::write(stream,req);

        beast::flat_buffer buffer;
        http::response<http::string_body> response;
        http::read(stream,buffer, response);
        stream.socket().shutdown(tcp::socket::shutdown_both);
        return json::parse(response.body());
    }

    void store_content(const std::string& server_address, const std::string& port,const std::string& path, std::vector<char> data){
        net::io_context ioc;
        auto stream = connect_stream(ioc,server_address,port);

        http::request<http::vector_body<char>> req{http::verb::post,server_address,http_version};
        req.set(http::field::host,server_address);
        req.set(http::field::user_agent, BOOST_BEAST_VERSION_STRING);
        req.body() = std::move(data);
        req.prepare_payload();
        http::write(stream,std::move(req));







    }




class RESTStorageClient : public StorageSpace::StreamProvider {
public:
    RESTStorageClient(const std::string& server_address, const std::string& group, const ISMRMRD::IsmrmrdHeader& header) : server_address(server_address), group(group),header(header) {
    }

    ~RESTStorageClient() override = default;

    std::vector<std::string> content(const std::string &key) const override {



    }

    std::vector<char> fetch(const std::string &key) const override {
        return std::unique_ptr<std::istream>();
    }
    void store(const std::string &key, const std::vector<char>& values) const override {
        return std::unique_ptr<std::ostream>();
    }

private:
    const std::string server_address;
    const std::string group;
    const ISMRMRD::IsmrmrdHeader& header
};
}