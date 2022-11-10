#include "StorageSetup.h"

#include <iterator>

#include <curl/curl.h>
#include <date/date.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using namespace Gadgetron::Storage;

namespace {
// Ensures that curl_global_init and curl_global_cleanup are called once
class CurlGlobalStateGuard {
  public:
    CurlGlobalStateGuard() { curl_global_init(CURL_GLOBAL_DEFAULT); }
    ~CurlGlobalStateGuard() { curl_global_cleanup(); }
};
static CurlGlobalStateGuard curl_global_state;

// Copies the response body
size_t content_write_callback(char* ptr, size_t size, size_t nmemb, std::stringstream* data) {
    data->write(ptr, size * nmemb);
    return size * nmemb;
}

// Copies the request body
size_t content_read_callback(char* dest, size_t size, size_t nmemb, std::istream* stream) {
    return stream->readsome(dest, size * nmemb);
}

template <typename T> using Handle = std::unique_ptr<T, std::function<void(T*)>>;

Handle<CURL> create_curl_handle() {
    Handle<CURL> handle(curl_easy_init(), curl_easy_cleanup);
    if (!handle) {
        throw std::runtime_error("unable to create CURL instance");
    }

    curl_easy_setopt(handle.get(), CURLOPT_NOPROGRESS, 1L);
    curl_easy_setopt(handle.get(), CURLOPT_MAXREDIRS, 50L);
    curl_easy_setopt(handle.get(), CURLOPT_TCP_KEEPALIVE, 1L);

    return handle;
}

StorageItem storage_item_from_json(json j) {
    StorageItem s;
    StorageItemTags::Builder tag_builder(j["subject"].get<std::string>());
    for (const auto& p : j.items()) {
        if (p.key() == "subject") {
            // Skip, we used that in the constructor
        } else if (p.key() == "device") {
            tag_builder.with_device(p.value().get<std::string>());
        } else if (p.key() == "session") {
            tag_builder.with_session(p.value().get<std::string>());
        } else if (p.key() == "name") {
            tag_builder.with_name(p.value().get<std::string>());
        } else if (p.key() == "contentType") {
            s.contentType = p.value().get<std::string>();
        } else if (p.key() == "data") {
            s.data = p.value().get<std::string>();
        } else if (p.key() == "location") {
            s.location = p.value().get<std::string>();
        } else if (p.key() == "lastModified") {
            std::istringstream in(p.value().get<std::string>());
            in >> date::parse("%FT%TZ", s.lastModified);
        } else if (p.key() == "expires") {
            std::istringstream in(p.value().get<std::string>());
            decltype(s.expires)::value_type expiration;
            in >> date::parse("%FT%TZ", expiration);
            s.expires = expiration;
        } else {
            if (p.value().is_array()) {
                for (auto& v : p.value()) {
                    tag_builder.with_custom_tag(p.key(), v.get<std::string>());
                }
            } else {
                tag_builder.with_custom_tag(p.key(), p.value().get<std::string>());
            }
        }
    }
    s.tags = tag_builder.build();
    return s;
}

std::string tags_to_query_string(StorageItemTags const& tags) {
    std::stringstream ss;
    ss << "?subject=" << tags.subject;
    if (tags.device.has_value()) {
        ss << "&device=" << *tags.device;
    }
    if (tags.session.has_value()) {
        ss << "&session=" << *tags.session;
    }
    if (tags.name.has_value()) {
        ss << "&name=" << *tags.name;
    }
    for (const auto& t : tags.custom_tags) {
        ss << "&" << t.first << "=" << t.second;
    }
    return ss.str();
}

StorageItemList get_item_list(std::string const& url) {
    std::stringstream response_body;
    auto curl_handle = create_curl_handle();
    curl_easy_setopt(curl_handle.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEFUNCTION, content_write_callback);
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEDATA, &response_body);

    auto res = curl_easy_perform(curl_handle.get());
    if (res != CURLE_OK) {
        throw std::runtime_error("Failed to list items: " + std::string(curl_easy_strerror(res)));
    }

    long status_code;
    curl_easy_getinfo(curl_handle.get(), CURLINFO_RESPONSE_CODE, &status_code);
    if (status_code != 200) {
        throw std::runtime_error("Storage server error when listing items.\n"
                                 "HTTP status code: " +
                                 std::to_string(status_code) + "\n" + "Body: " + response_body.str());
    }

    json j = json::parse(response_body.str());
    StorageItemList list;
    for (auto item : j["items"]) {
        list.items.push_back(storage_item_from_json(item));
    }

    if (j.contains("nextLink")) {
        list.complete = false;
        list.continuation = j["nextLink"];
    } else {
        list.complete = true;
    }

    return list;
}

} // namespace

namespace Gadgetron::Storage {

StorageItemList StorageClient::list_items(StorageItemTags const& tags, size_t limit) {
    std::string query_string = tags_to_query_string(tags);
    query_string += "&_limit=" + std::to_string(limit);

    std::string url = base_url + "/v1/blobs" + query_string;
    return get_item_list(url);
}

StorageItemList StorageClient::get_next_page_of_items(StorageItemList const& page) {
    if (page.complete || page.continuation.empty()) {
        return StorageItemList{ .complete = true };
    }

    return get_item_list(page.continuation);
}

std::shared_ptr<std::istream> StorageClient::get_latest_item(StorageItemTags const& tags) {
    std::string url = base_url + "/v1/blobs/data/latest" + tags_to_query_string(tags);
    return get_item_by_url(url);
}

std::shared_ptr<std::istream> StorageClient::get_item_by_url(const std::string& url) {
    // Note that we are loading the entire response into memory here, which will be 
    // problematic for large items.

    auto response_body = std::make_shared<std::stringstream>();
    auto curl_handle = create_curl_handle();
    curl_easy_setopt(curl_handle.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEFUNCTION, content_write_callback);
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEDATA, response_body.get());

    auto res = curl_easy_perform(curl_handle.get());
    if (res != CURLE_OK) {
        throw std::runtime_error("Failed to get item: " + std::string(curl_easy_strerror(res)));
    }

    long status_code;
    curl_easy_getinfo(curl_handle.get(), CURLINFO_RESPONSE_CODE, &status_code);
    if (status_code == 404) {
        return {};
    }

    if (status_code != 200) {
        throw std::runtime_error("Storage server error when getting item.\n"
                                 "HTTP status code: " +
                                 std::to_string(status_code) + "\n" + "Body: " + response_body->str());
    }

    return response_body;
}

StorageItem StorageClient::store_item(StorageItemTags const& tags, std::istream& data,
                                      std::optional<std::chrono::seconds> time_to_live) {
    auto query_string = tags_to_query_string(tags);
    if (time_to_live) {
        query_string += "&_ttl=" + std::to_string(time_to_live->count()) + "s";
    }

    std::string url = base_url + "/v1/blobs/data" + query_string;
    std::stringstream response_body;

    auto curl_handle = create_curl_handle();
    curl_easy_setopt(curl_handle.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl_handle.get(), CURLOPT_POST, 1L);
    curl_easy_setopt(curl_handle.get(), CURLOPT_READFUNCTION, content_read_callback);
    curl_easy_setopt(curl_handle.get(), CURLOPT_READDATA, &data);
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEFUNCTION, content_write_callback);
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEDATA, &response_body);

    // use chunked transfer encoding and disable Expect: 100-continue
    Handle<curl_slist> headers_handle(
        curl_slist_append(curl_slist_append(nullptr, "Transfer-Encoding: chunked"), "Expect:"), curl_slist_free_all);
    curl_easy_setopt(curl_handle.get(), CURLOPT_HTTPHEADER, headers_handle.get());

    auto res = curl_easy_perform(curl_handle.get());
    if (res != CURLE_OK) {
        throw std::runtime_error("Failed to store item: " + std::string(curl_easy_strerror(res)));
    }

    long status_code;
    curl_easy_getinfo(curl_handle.get(), CURLINFO_RESPONSE_CODE, &status_code);
    if (status_code != 201) {
        throw std::runtime_error("Storage server error when storing item.\n"
                                 "HTTP status code: " +
                                 std::to_string(status_code) + "\n" + "Body: " + response_body.str());
    }

    return storage_item_from_json(json::parse(response_body.str()));
}

std::optional<std::string> StorageClient::health_check() {
    std::string url = base_url + "/healthcheck";
    std::stringstream response_body;

    auto curl_handle = create_curl_handle();
    curl_easy_setopt(curl_handle.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEFUNCTION, content_write_callback);
    curl_easy_setopt(curl_handle.get(), CURLOPT_WRITEDATA, &response_body);

    auto res = curl_easy_perform(curl_handle.get());
    if (res != CURLE_OK) {
        return "Failed to perform health check: " + std::string(curl_easy_strerror(res));
    }

    long status_code;
    curl_easy_getinfo(curl_handle.get(), CURLINFO_RESPONSE_CODE, &status_code);
    if (status_code == 200) {
        return std::nullopt;
    }

    return "Storage server error when performing health check.\n"
           "HTTP status code: " +
           std::to_string(status_code) + "\n" + "Body: " + response_body.str();
}

} // namespace Gadgetron::Storage
