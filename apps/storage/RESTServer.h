#pragma once
#include <boost/beast.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/spawn.hpp>
#include <boost/hana.hpp>
#include <iostream>
#include <future>

namespace Gadgetron::Storage::REST {
    namespace beast = boost::beast;         // from <boost/beast.hpp>
    namespace http = beast::http;           // from <boost/beast/http.hpp>
    namespace net = boost::asio;            // from <boost/asio.hpp>
    using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>
    namespace asio = boost::asio;
    namespace hana = boost::hana;


    auto make_failure_response(std::string_view message, http::status status){
        http::response<http::string_body> res{status, 10 };
        res.set(http::field::server, BOOST_BEAST_VERSION_STRING);
        res.set(http::field::content_type, "text/plain");
        res.keep_alive(false);
        res.body() = message;
        res.prepare_payload();
        return res;
    }

    template<class... Endpoints>

    class RESTEndpoints {
    public:
        RESTEndpoints(Endpoints &&... endpoints) : endpoints(std::forward<Endpoints>(endpoints)...) {}

        template<class SyncReadStream, class Send, class Buffer >
        beast::error_code handle_request(SyncReadStream &stream,Send &&send, Buffer& buffer, asio::yield_context yield ) {
            beast::error_code ec;
            http::request_parser<http::empty_body> req_parser;
            req_parser.body_limit(128ull * 1024ull * 1024ull * 1024ull); //We support files up to 128GB. For now.
            http::async_read_header(stream, buffer, req_parser,yield[ec]);
            if (ec == http::error::end_of_stream) return ec;
            if (ec) return send(make_failure_response(ec.message(),http::status::bad_request));

            auto req = req_parser.get();

            auto reader = [&stream, &buffer, yield](auto&& request_parser){
                beast::error_code  errorCode;
                http::async_read(stream,buffer,request_parser,yield[errorCode]);
                return errorCode;
            };

            std::string fail_message = "The resource '" + std::string(req.target()) + "' was not found.";

            try {
                bool accepted = hana::fold_left(endpoints, false, [&](bool state, auto &&endpoint) {
                    if (state) return true;
                    bool accepts = endpoint.accept(req);
                    if (!accepts) return false;
                    ec = endpoint.handle(reader, send, std::move(req_parser));
                    return true;
                });

                if (accepted) return ec;
            } catch (const std::exception& error){
                fail_message = error.what();
            }

            ec = send(make_failure_response(fail_message,http::status::bad_request));
            return ec;

        }

    private:

        hana::tuple<Endpoints...> endpoints;

    };

namespace internal{

    template< bool isRequest, class Body, class Fields>
    auto create_serializer(http::message<isRequest,Body,Fields> & message){
        return http::serializer<isRequest,Body,Fields>(message);
    }
}

template<class RequestHandler>
void handle_session( beast::tcp_stream stream, asio::yield_context yield, RequestHandler& handler ) {
    try {
        bool close = false;
        beast::flat_buffer buffer;
        auto sender = [&stream, &close, yield](auto &&message) {
            beast::error_code ec;
            close = message.need_eof();
            auto sr = internal::create_serializer(message);
            http::async_write(stream, sr, yield[ec]);
            return ec;
        };

        for (;;) {
            beast::error_code ec = handler.handle_request(stream, sender, buffer, yield);
            if (ec == http::error::end_of_stream) break;
            if (close) break;
        }
        stream.socket().shutdown(tcp::socket::shutdown_send);
    } catch(const std::runtime_error& error){
        std::cerr <<"HTTP session failed with error " <<  error.what() << std::endl;
    }

}


template<class RequestHandler>
void navi( asio::io_context& ioc, tcp::endpoint endpoint,  RequestHandler& handler,std::promise<unsigned short>& port, asio::yield_context yield){
    tcp::acceptor acceptor(ioc);
    acceptor.open(endpoint.protocol());
    acceptor.set_option(asio::socket_base::reuse_address(true));
    acceptor.bind(endpoint);
    acceptor.listen(asio::socket_base::max_listen_connections);
    port.set_value(acceptor.local_endpoint().port());
    while (!ioc.stopped()){
        tcp::socket socket(ioc);
        acceptor.async_accept(socket,yield);
        asio::spawn(acceptor.get_executor(),[&socket,&handler](auto yield){ handle_session(beast::tcp_stream(std::move(socket)),yield,handler);});
    }
}



}

