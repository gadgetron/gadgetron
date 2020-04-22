//
// Created by david on 17/04/2020.
//

#include "storage_server.h"
#include "Context.h"
#include "Storage.h"
#include <boost/process.hpp>
#include <cpprest/http_client.h>
#include <cpprest/interopstream.h>

Gadgetron::Server::StorageServer Gadgetron::Server::start_storage_server(
    const boost::filesystem::path& working_directory) {
    namespace bp = boost::process;

    auto process = bp::child(
        bp::search_path("python3"), "storage_server.py", "--storage_dir", working_directory.string(), "--port", "9102");
    return { std::move(process), "localhost:9102" };
}

namespace {
    using namespace Gadgetron::Core;
    using namespace utility;           // Common utilities like string conversions
    using namespace web;               // Common features like URIs.
    using namespace web::http;         // Common HTTP functionality
    using namespace web::http::client; // HTTP client features
    using namespace concurrency::streams;

    class RestStorage : public StorageSpace::StreamProvider {
    public:
        RestStorage(const std::string& server_address, const std::string& group)
            : server_address(server_address), group(group) { }

        std::unique_ptr<std::istream> fetch(const std::string& key) const override {
            http_client client(U(server_address));
            uri_builder builder;
            builder.append(U("/v1/"));
            builder.append(U(group));
            builder.append(U(key));
            auto data = client.request(methods::GET, builder.to_string())
                            .then([=](http_response response) { return response.extract_json(); })
                            .then([=](json::value node) mutable {
                                auto blob_uri = node["contents"][0]["uri"].as_string();
                                return client.request(methods::GET, blob_uri);
                            })
                            .then([=](http_response response) { return new async_istream<char>(response.body()); });
            return std::unique_ptr<std::istream>(data.get());
        }

        bool contains(const std::string& key) const override {
            http_client client(U(server_address));
            uri_builder builder;
            builder.append(U("/v1/"));
            builder.append(U(group));
            builder.append(U(key));
            auto json_response = client.request(methods::GET, builder.to_string())
                                     .then([=](http_response response) {
                                         return response.extract_json();
                                     })
                                     .get();

            return json_response["contents"].size() > 0;
        }

        void store(const std::string& key, std::istream& data) const override {

            http_client client(U(server_address));
            uri_builder builder;
            builder.append(U("/v1/blobs"));

            json::value message  = json::value::object();
            message["operation"] = json::value::string("push");
            message["arguments"] = json::value::array(std::vector<json::value>{ json::value::string(key) });

            auto response = client.request(methods::PATCH, builder.to_string(), message)
                                .then([=](http_response response) { return response.extract_json(); })
                                .then([=](json::value node) {
                                    auto id = node["contents"][0]["id"].as_string();
                                    uri_builder builder;
                                    builder.append("/blobs/");
                                    builder.append(id);
                                    return builder.to_string();
                                });
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
