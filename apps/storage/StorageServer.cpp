//
// Created by dch on 4/9/21.
//

#include "StorageServer.h"
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
#include <range/v3/algorithm.hpp>


using namespace Gadgetron::Storage;
namespace beast = boost::beast;         // from <boost/beast.hpp>
namespace http = beast::http;           // from <boost/beast/http.hpp>
namespace net = boost::asio;            // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>
namespace asio = boost::asio;

using json = nlohmann::json;

using namespace Gadgetron::Storage::DB;


namespace {
    enum class StorageSpace {
        session,
        scanner,
        measurement
    };

    static std::vector<std::pair<std::string, StorageSpace>> storage_space_names = {{"session", StorageSpace::session},
                                                                                    {"scanner", StorageSpace::scanner},
                                                                                    {"measurement", StorageSpace::measurement}};

    void from_json(const json &j, StorageSpace &storagespace) {
        auto string_val = j.get<std::string>();
        auto it = std::find_if(storage_space_names.begin(), storage_space_names.end(),
                               [&](auto &val) { return val.first == string_val; });
        if (it == storage_space_names.end()) throw std::runtime_error("Invalid storage space provided");
        storagespace = it->second;
    }

    void to_json(json &j, const StorageSpace &storagespace) {
        auto it = std::find_if(storage_space_names.begin(), storage_space_names.end(),
                               [&](auto &val) { return val.second == storagespace; });
        if (it == storage_space_names.end()) throw std::runtime_error("Infernal server error");
        j = it->first;
    }


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
        boost::filesystem::path blob_folder;

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
            req2_parser.body_limit(128ull * 1024ull * 1024ull * 1024ull); //We support files up to 128GB. For now.
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
        boost::filesystem::path blob_folder;


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
            BOOST_HANA_DEFINE_STRUCT(RetrievalResponse, (std::string, storage_path),
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
            auto key = full_key(retrieval_request.storagespace, retrieval_request.subject, retrieval_request.key);
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
            BOOST_HANA_DEFINE_STRUCT(WriteResponse, (std::string, storage_path));
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

                auto pending_write = PendingWrite{key, deadline, {blob_id, now, now + write_request.storage_duration}};

                db->pending_writes.set(to_string(blob_id), pending_write);
                return send(json_response(make_response(blob_id), req2.get()));
            } catch (const std::exception &error) {
                std::string error_message = error.what();
                return send(string_response(error_message, http::status::bad_request, req2.get()));
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

void cleanup(const boost::system::error_code &ec, asio::deadline_timer &timer, std::shared_ptr<DB::Database> db, boost::filesystem::path blob_folder) {

    auto now = boost::posix_time::second_clock::universal_time();

    using namespace ranges;
    auto stale_writes = db->pending_writes | views::filter([&](const auto& key_value) {
        const auto&[key, pending] = key_value;
        return pending.transaction_expiration < now;
    }) | views::keys;

    db->pending_writes.delete_keys(stale_writes);

    auto is_blob_stale = [&](const auto& meta ){
        return meta.deletion_time < now;
    };

    std::vector<std::pair<std::string, std::vector<BlobMeta>>> updated_lists;
    std::vector<boost::uuids::uuid> blobs_to_delete;
    for (auto [key,meta_list] : db->blobs ){
        if (!any_of(meta_list,is_blob_stale)) continue;
        auto new_list = meta_list | views::filter([&](const auto& meta ){return !is_blob_stale(meta);}) | to<std::vector>();
        updated_lists.emplace_back(key,std::move(new_list));

        auto marked_for_deletion = meta_list | views::filter(is_blob_stale) | views::transform([](const BlobMeta& meta){
            return meta.blob_id;
        });
        blobs_to_delete.insert(blobs_to_delete.end(),begin(marked_for_deletion),end(marked_for_deletion));
    }


    db->blobs.update(updated_lists);

    //Give people 30 minutes grace period to actually retrieve the data.
    for (const auto& blob_id : blobs_to_delete){
        db->scheduled_deletions.set(to_string(blob_id),{blob_id,now+boost::posix_time::time_duration(0,30,0)});
    }

    for (const auto& [id,deletion_record] : db->scheduled_deletions ){
        if (deletion_record.deletion_time > now) continue;
        auto path = blob_folder / to_string(deletion_record.blob_id);
        boost::system::error_code ec;
        boost::filesystem::remove(path,ec);
        db->scheduled_deletions.delete_key(id);
    }

    timer.expires_at(timer.expires_at() + boost::posix_time::minutes(5));
    timer.async_wait([&timer, db,blob_folder](const auto &ec) { cleanup(ec, timer, db,std::move(blob_folder)); });
}

static void ensure_exists(const boost::filesystem::path &folder) {
    if (exists(folder)) {
        if (!is_directory(folder)) throw std::runtime_error("Specified path " + folder.string() + " must be a folder");
        return;
    }

    create_directories(folder);
}

Gadgetron::Storage::StorageServer::StorageServer(unsigned short port, const boost::filesystem::path &database_folder,
                                                 const boost::filesystem::path &blob_folder) : ioContext{} {

    ensure_exists(database_folder);
    ensure_exists(blob_folder);
    std::promise<unsigned short> bound_port_promise;

    this->server_thread = std::thread([this, database_folder, blob_folder, port, &bound_port_promise] {

        auto database = DB::Database::make_db(database_folder);
        REST::RESTEndpoints endpoints(InfoEndPoint{}, BlobStorageEndPoint{database, blob_folder},
                                      BlobRetrievalEndPoint{database, blob_folder}, WriteRequestEndPoint{database},
                                      DataInfoEndPoint{database});

        asio::spawn(ioContext, [this, &endpoints, port, &bound_port_promise](auto yield) {
            REST::navi(ioContext, tcp::endpoint(asio::ip::tcp::v6(), port), endpoints, bound_port_promise, yield);
        });

        asio::deadline_timer timer(ioContext, boost::posix_time::seconds(1));
        timer.async_wait([&](const auto& ec){ cleanup(ec,timer,database,blob_folder);});
        ioContext.run();
    });

    this->bound_port = bound_port_promise.get_future().get();
}

StorageServer::~StorageServer() {
    ioContext.stop();
    if (server_thread.joinable())
        server_thread.join();

}

unsigned short StorageServer::port() {
    return this->bound_port;
}

void StorageServer::run_forever() {
    this->server_thread.join();
}
