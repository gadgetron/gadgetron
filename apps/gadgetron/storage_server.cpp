//
// Created by david on 17/04/2020.
//

#include "storage_server.h"
#include "Context.h"
#include "Storage.h"
#include <boost/process.hpp>
#include <cpprest/http_client.h>
#include <cpprest/interopstream.h>
#include <cpprest/producerconsumerstream.h>
#include "gadgetron_config.h"

Gadgetron::Server::StorageServer Gadgetron::Server::start_storage_server(
    const boost::filesystem::path& working_directory) {
    namespace bp = boost::process;

    auto process = bp::child(
        GADGETRON_PYTHON_EXECUTABLE, "storage_server.py", "--storage-dir", working_directory.string(), "--port", "9102");
    return { std::move(process), "http://localhost:9102" };
}

namespace {
    using namespace Gadgetron::Core;
    using namespace utility;           // Common utilities like string conversions
	using namespace utility::conversions;
    using namespace web;               // Common features like URIs.
    using namespace web::http;         // Common HTTP functionality
    using namespace web::http::client; // HTTP client features
    using namespace concurrency::streams;

    class storage_iostream : public async_iostream<char> {
        public:

        storage_iostream(const std::string& key, const std::string& group, const std::string& server_address, producer_consumer_buffer<char> buffer = producer_consumer_buffer<char>{}) : async_iostream<char>(buffer),  internal_buffer(buffer){

            http_client restclient(to_string_t(server_address));
            uri_builder builder;
            builder.append(U("/v1/blobs"));

            auto other_buffer = streambuf<uint8_t>(internal_buffer);

			

            task = restclient.request(methods::PUT,builder.to_string(),other_buffer.create_istream())
                .then([=](http_response response){ return response.extract_json();})
                .then([=](json::value node){
                    auto id = node[U("id")].as_string();
                    return id;
                })
                .then([=](string_t id) mutable {
                    json::value message  = json::value::object();
                    message[U("operation")] = json::value::string(U("push"));
                    message[U("arguments")] = json::value::array(std::vector<json::value>{ json::value::string(to_string_t(id)) });
                    uri_builder builder(U("/v1"));
                    builder.append(to_string_t(group));
                    builder.append(to_string_t(key));

                    return restclient.request(methods::PATCH,builder.to_string(),message);
                })
                .then([=](http_response response){
                    return response.extract_json();
                })
                .then([=](json::value node){
					
                });

        }

        ~storage_iostream(){
            internal_buffer.sync().wait();
            internal_buffer.close(std::ios_base::out).wait();
            
            task.wait();

        }

        private: 

        producer_consumer_buffer<char> internal_buffer;
        
        pplx::task<void> task;


    };
   



    class RestStorage : public StorageSpace::StreamProvider {
    public:
        RestStorage(const std::string& server_address, const std::string& group)
            : server_address(server_address), group(group) { }

        std::unique_ptr<std::istream> fetch(const std::string& key) const override {
            http_client client(to_string_t(server_address));
            uri_builder builder;
            builder.append(U("/v1/"));
            builder.append(to_string_t(group));
            builder.append(to_string_t(key));
            auto data = client.request(methods::GET, builder.to_string())
                            .then([=](http_response response) { return response.extract_json(); })
                            .then([=](json::value node) mutable {
                                auto blob_uri = node[U("contents")][0][U("uri")].as_string();
                                return client.request(methods::GET, blob_uri);
                            })
                            .then([=](http_response response) { return new async_istream<char>(response.body()); });
            return std::unique_ptr<std::istream>(data.get());
        }

        bool contains(const std::string& key) const override {
            http_client client(to_string_t(server_address));
            uri_builder builder;
            builder.append(U("/v1/"));
            builder.append(to_string_t(group));
            builder.append(to_string_t(key));
            auto json_response = client.request(methods::GET, builder.to_string())
                                     .then([=](http_response response) {
                                         return response.extract_json();
                                     })
                                     .get();


            return json_response.size() > 0;
        }

        std::unique_ptr<std::ostream> store(const std::string& key) const override {
            return std::make_unique<storage_iostream>(key,group,server_address);
        }

    private:
        std::string server_address;
        std::string group;
    };
}

Gadgetron::Core::Storage Gadgetron::Server::setup_storage(
    const StorageServer::Address& address, const Core::Context::Header& header) {

    using namespace std::string_literals;

    std::string session_uid;
    if (header.measurementInformation.is_present()) {
        session_uid = header.measurementInformation->seriesInstanceUIDRoot.is_present()
                          ? header.measurementInformation->seriesInstanceUIDRoot.get()
                          : "default";
    } else {
        session_uid = "default";
    }

  return Core::Storage{
      Core::StorageSpace(std::make_shared<RestStorage>(address,"session/"s + session_uid)),
      Core::StorageSpace(std::make_shared<RestStorage>(address,"noise")),
      Core::StorageSpace(std::make_shared<RestStorage>(address,"debug")),
  };
}
