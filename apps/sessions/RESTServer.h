#pragma once
#include <boost/beast.hpp>
#include <boost/hana.hpp>


namespace Gadgetron::Sessions::REST {
    namespace beast = boost::beast;         // from <boost/beast.hpp>
    namespace http = beast::http;           // from <boost/beast/http.hpp>
    namespace net = boost::asio;            // from <boost/asio.hpp>
    using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>

    using bh = boost::hana;


    template<class... Endpoints>

    class RESTServer {

        RESTServer(Endpoints...&& endpoints) : endpoints(std::forward<Endpoints>(endpoints)...){}

        template<class SyncReadStream, class Send>
       void handle_request(SyncReadStream& stream, Send&& send){
            http::request_parser<empty_body> req;
            http::read_header(stream, req);

            hana::fold_left(endpoints,false, [&](auto&& endpoint, bool state){
                if (state) return true;
                bool accepts = endpoint.accept(req);
                if (!accepts) return false;
                send(endpoint.handle(stream,std::move(req)));
                return true;
            } );



        }




    private:

        hana::tuple<Endpoints> endpoints;

    };
}


