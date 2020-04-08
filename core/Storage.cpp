#include "Storage.h"
#include <cpprest/http_client.h>
#include <cpprest/interopstream.h>

namespace Gadgetron::Core {
    StorageSpace::StorageSpace(std::unique_ptr<StreamProvider> provider) : provider(std::move(provider)) {}
}


namespace {
    using namespace Gadgetron::Core;
    using namespace utility;                    // Common utilities like string conversions
    using namespace web;                        // Common features like URIs.
    using namespace web::http;                  // Common HTTP functionality
    using namespace web::http::client;          // HTTP client features
    using namespace concurrency::streams;



    class RestStorage : public StorageSpace::StreamProvider {
    public:
        RestStorage(const std::string& server_address) : client(U(server_address)){

        }

        std::unique_ptr<std::istream> fetch(const std::string &key) override {
            
            uri_builder builder;
            builder.append(U("/sessions/"));
            builder.append(U(key));
            auto data = client.request(methods::GET,builder.to_string())
                    .then([=](http_response response){
                        return response.extract_json();
                    })
                    .then([=](json::value node){
                        auto id = node["contents"][0]["id"].as_string();
                        uri_builder builder;
                        builder.append("/blobs/");
                        builder.append(id);
                        return builder.to_string();
                    })
                    .then([=](const std::string& blob_uri){
                       return client.request(methods::GET,blob_uri);
                    })
                    .then([=](http_response response){
                        return std::make_unique<async_istream<char>>(response.body());
                    })
                    ;
            return data.get();

        }

        void store(const std::string &key, const std::function<void(std::ostream&)>& data) override {
            uri_builder builder;
            builder.append(U("/sessions"));
            builder.append(U("/blobs"));
            ;

            client.request(methods::POST,builder.to_string(),)
            
            json::value message = json::value::object();
            message["operation"]= json::value::string("push");
            message["arguments"]= json::value::array(std::vector<json::value>{json::value::string(key)});

            auto response = client.request(methods::PATCH,builder.to_string(),message)
                    .then([=](http_response response){
                        return response.extract_json();
                    })
                    .then([=](json::value node )){
                        auto id = node["contents"][0]["id"].as_string();
                        uri_builder builder;
                        builder.append("/blobs/");
                        builder.append(id);
                        return builder.to_string();

                    }
            ;


        }

    private:
        http_client client;


    };
}