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
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/string_generator.hpp>
#include "Database.h"

using namespace Gadgetron::Sessions;
namespace beast = boost::beast;         // from <boost/beast.hpp>
namespace http = beast::http;           // from <boost/beast/http.hpp>
namespace net = boost::asio;            // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>
namespace asio = boost::asio;

using json = nlohmann::json;

using namespace Gadgetron::Sessions::DB;


namespace {
    enum class StorageSpace {
        session,
        scanner,
        debug
    };

    NLOHMANN_JSON_SERIALIZE_ENUM(StorageSpace, {
        { StorageSpace::session, "session" },
        { StorageSpace::scanner, "scanner" },
        { StorageSpace::debug, "debug" }
    }
    );

    std::string name(StorageSpace space) {
        return json(space).get<std::string>();
    }


    std::string full_key(StorageSpace space, std::string_view subject, std::string_view key) {

        auto escape_string = [](const auto &str) {
            std::string result;
            if (str.size() > std::numeric_limits<uint32_t>::max()) throw std::runtime_error("Key length out of bounds");

            uint32_t length = str.size();

            result.reserve(sizeof(uint32_t) + str.size());
            result.append(reinterpret_cast<const char *>(&length), sizeof(length));
            result.append(str);
            return result;
        };


        std::string result = name(space) + '/';
        result.append(escape_string(subject));
        result.push_back('/');
        result.append(key);
        return result;
    }


    template<class REQ>
    auto json_response(const json &object, const REQ &req) {
        http::response<http::string_body> response(http::status::ok, req.version());
        response.set(http::field::server, BOOST_BEAST_VERSION_STRING);
        response.set(http::field::content_type, "application/json");
        response.keep_alive(req.keep_alive());
        response.body() = object.dump();
        response.prepare_payload();
        return response;
    }

    template<class REQ>
    auto string_response(std::string_view message, http::status status, const REQ &req) {
        http::response<http::string_body> response(status, req.version());
        response.set(http::field::server, BOOST_BEAST_VERSION_STRING);
        response.set(http::field::content_type, "text/plain");
        response.keep_alive(req.keep_alive());
        response.body() = message;
        response.prepare_payload();
        return response;
    }

    template<class REQ>
    auto empty_response(http::status status, const REQ &req) {
        http::response<http::string_body> response(status, req.version());
        response.set(http::field::server, BOOST_BEAST_VERSION_STRING);
        response.keep_alive(req.keep_alive());
        response.prepare_payload();
        return response;
    }

    struct InfoEndPoint {

        bool accept(const http::request<http::empty_body> &req) {
            if (req.method() == http::verb::get && req.target() == "/v1/info") return true;
            return false;
        }

        template<class Read, class Send>
        beast::error_code handle(Read &read, Send &&send, http::request_parser<http::empty_body> &&req) {
            json version = {{"version", "0.01alpha"}};
            return send(json_response(version, req.get()));
        }
    };

    template<class REQ>
    std::string get_blob_key(const REQ &req, std::string_view endpoint) {
        auto target = req.target();
        target.remove_prefix(endpoint.size());
        return std::string(target);
    }

    struct BlobStorageEndPoint {

        std::shared_ptr<Database> database;
        std::filesystem::path blob_folder;

        bool accept(const http::request<http::empty_body> &req) {
            if (req.method() == http::verb::patch && req.target().starts_with(endpoint)) return true;
            return false;
        }

        template<class Read, class Send>
        beast::error_code handle(Read &read, Send &&send, http::request_parser<http::empty_body> &&req) {
            beast::error_code ec;

            auto key = get_blob_key(req.get(), endpoint);

            auto potential_meta = database->pending_writes[key];

            if (!potential_meta) {
                return send(string_response("Unknown ID", http::status::bad_request, req.get()));
            }

            auto &meta = *potential_meta;

            auto path = blob_folder / to_string(meta.meta.blob_id);
            http::request_parser<http::file_body> req2_parser(std::move(req));
            req2_parser.body_limit(128ull*1024ull*1024ull*1024ull); //We support files up to 128GB. For now.
            req2_parser.get().body().open(path.c_str(), beast::file_mode::write, ec);
            if (ec) {
                return send(string_response(ec.message(), http::status::insufficient_storage, req2_parser.get()));
            }

            auto ec2 = read(req2_parser);
            if (ec2) return ec2;


            database->pending_writes.delete_key(key);
            database->blobs.push_back(meta.key, meta.meta);

            return send(empty_response(http::status::ok, req2_parser.get()));
        }

        static constexpr const char *endpoint = "/v1/blobs/";

    };

    struct BlobRetrievalEndPoint {

        std::shared_ptr<Database> database;
        std::filesystem::path blob_folder;


        bool accept(const http::request<http::empty_body> &req) {
            if (req.method() == http::verb::get && req.target().starts_with(endpoint)) return true;
            return false;
        }

        template<class Read, class Send>
        beast::error_code handle(Read &read, Send &&send, http::request_parser<http::empty_body> &&req) {

            boost::uuids::string_generator gen;
            boost::uuids::uuid blob_id = gen(get_blob_key(req.get(), endpoint));

            auto path = blob_folder / to_string(blob_id);

            http::response<http::file_body> response(http::status::ok, req.get().version());
            response.set(http::field::server, BOOST_BEAST_VERSION_STRING);
            response.set(http::field::content_type, "application/json");
            beast::error_code ec;
            response.body().open(path.c_str(), beast::file_mode::read, ec);
            if (ec) {
                return send(string_response(ec.message(), http::status::not_found, req.get()));
            }

            response.prepare_payload();
            return send(response);
        }

        static constexpr const char *endpoint = BlobStorageEndPoint::endpoint;

    };

    struct DataInfoEndPoint {

        struct RetrievalRequest {
            BOOST_HANA_DEFINE_STRUCT(RetrievalRequest, (StorageSpace, storagespace), (std::string, subject),
                                     (std::string, key));
        };

        struct RetrievalResponse {
            BOOST_HANA_DEFINE_STRUCT(RetrievalResponse, (std::string, storagepath),
                                     (boost::posix_time::ptime, creation_time),
                                     (boost::posix_time::ptime, deletion_time));
        };

        bool accept(const http::request<http::empty_body> &req) {
            if (req.method() == http::verb::get && req.target().starts_with(endpoint)) return true;
            return false;
        }

        template<class Read, class Send>
        beast::error_code handle(Read &read, Send &&send, http::request_parser<http::empty_body> &&req) {
            http::request_parser<http::string_body> req2(std::move(req));
            if (auto ec = read(req2)) return ec;
            auto retrieval_request = json::parse(req2.get().body()).get<RetrievalRequest>();
            auto key = full_key(retrieval_request.storagespace,retrieval_request.subject, retrieval_request.key);
            std::cout << "Retrieving from key " << key << std::endl;
            auto blobs = database->blobs[key];

            auto response_blobs = ranges::reverse_view(blobs) | ranges::views::transform(
                    [](auto &&meta) -> RetrievalResponse {
                        return {std::string(BlobRetrievalEndPoint::endpoint) + to_string(meta.blob_id),
                                meta.creation_time, meta.deletion_time};
                    }) | ranges::to<std::vector>();

            return send(json_response(response_blobs, req2.get()));
        }


        static constexpr const char *endpoint = "/v1/data";
        std::shared_ptr<Database> database;
    };

    struct WriteRequestEndPoint {

        struct WriteRequest {
            BOOST_HANA_DEFINE_STRUCT(WriteRequest, (StorageSpace, storagespace), (std::string, subject),
                                     (std::string, key),
                                     (boost::posix_time::time_duration, storage_duration));
        };

        struct WriteResponse {
            BOOST_HANA_DEFINE_STRUCT(WriteResponse, (std::string, blob_path));
        };

        bool accept(const http::request<http::empty_body> &req) {
            if (req.method() == http::verb::post && req.target() == endpoint) return true;
            return false;
        }

        template<class Read, class Send>
        beast::error_code handle(Read &read, Send &&send, http::request_parser<http::empty_body> &&req) {

            http::request_parser<http::string_body> req2(std::move(req));
            auto ec = read(req2);

            if (ec) return ec;
            auto body = req2.get().body();
            try {
                auto j = json::parse(body);
                auto write_request = j.get<WriteRequest>();

                const std::string write_uuid = to_string(uuid_generator());
                const auto blob_id = uuid_generator();

                auto now = boost::posix_time::second_clock::universal_time();
                auto deadline = now + boost::posix_time::minutes(30);

                auto key = full_key(write_request.storagespace, write_request.subject, write_request.key);
                std::cout << "Storing key under: " << key << std::endl;

                auto pending_write = PendingWrite{key, deadline, {blob_id, now, now + write_request.storage_duration}};

                db->pending_writes.set(to_string(blob_id), pending_write);
                return send(json_response(make_response(blob_id), req2.get()));
            } catch (const std::exception& error){
                std::string error_message = error.what();
                return send(string_response(error_message,http::status::bad_request,req2.get()));
            }
        }

        std::shared_ptr<Database> db;

        boost::uuids::random_generator uuid_generator = {};
        static constexpr const char *endpoint = DataInfoEndPoint::endpoint;
    private:
        WriteResponse make_response(const boost::uuids::uuid &write_uuid) {

            return {std::string(BlobStorageEndPoint::endpoint) + to_string(write_uuid)};
        }

    };
}

void cleanup(const boost::system::error_code &ec, asio::deadline_timer &timer) {
    std::cout << "Cleaning up all da thing" << std::endl;

    timer.expires_at(timer.expires_at() + boost::posix_time::minutes(5));
    timer.async_wait([&](const auto &ec) { cleanup(ec, timer); });
}

static void ensure_exists(const std::filesystem::path& folder){
    if (exists(folder)){
        if (!is_directory(folder)) throw std::runtime_error("Specified path " + folder.string() + " must be a folder");
        return;
    }

    create_directories(folder);
}

Gadgetron::Sessions::SessionServer::SessionServer(unsigned short port, const std::filesystem::path &database_folder,
                                                  const std::filesystem::path &blob_folder): ioContext{} {

    ensure_exists(database_folder);
    ensure_exists(blob_folder);
    std::promise<unsigned short> bound_port_promise;

    this->server_thread = std::thread([this,database_folder, blob_folder, port,&bound_port_promise ] {

        auto database = std::make_shared<DB::Database>(database_folder);
        REST::RESTEndpoints endpoints(InfoEndPoint{}, BlobStorageEndPoint{database, blob_folder},
                                      BlobRetrievalEndPoint{database, blob_folder}, WriteRequestEndPoint{database},
                                      DataInfoEndPoint{database});

        asio::spawn(ioContext, [this, &endpoints, port, &bound_port_promise](auto yield) {
            REST::navi(ioContext, tcp::endpoint(asio::ip::tcp::v4(), port), endpoints,bound_port_promise, yield);
        });

        asio::deadline_timer timer(ioContext, boost::posix_time::seconds(1));
        //timer.async_wait([&](const auto& ec){ cleanup(ec,timer);});
        ioContext.run();
    });

    this->bound_port = bound_port_promise.get_future().get();
}

SessionServer::~SessionServer() {
    ioContext.stop();
    this->server_thread.join();

}

unsigned short SessionServer::port() {
    return this->bound_port;
}
