#include "Storage.h"

#include <iterator>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <cpr/cpr.h>
#include <date/date.h>
#include <nlohmann/json.hpp>

namespace bio = boost::iostreams;
using json = nlohmann::json;

namespace Gadgetron::Storage {

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

cpr::Parameters tags_to_query_parameters(StorageItemTags const& tags) {
    cpr::Parameters parameters{{"subject", tags.subject}};

    if (tags.device) {
        parameters.Add({"device", *tags.device});
    }
    if (tags.session) {
        parameters.Add({"session", *tags.session});
    }
    if (tags.name) {
        parameters.Add({"name", *tags.name});
    }
    for (auto const& [key, value] : tags.custom_tags) {
        parameters.Add({key, value});
    }

    return parameters;
}

std::string get_response_error_message(cpr::Response const& resp) {
    std::string details;
    if (resp.error.code != cpr::ErrorCode::OK) {
        if (resp.error.message.empty()) {
            return "CPR error code " + std::to_string((int)resp.error.code);
        } else {
            return resp.error.message;
        }
    }

    return resp.status_line + "\n" + resp.text;
}

StorageItemList read_items_from_response(cpr::Response const& resp) {
    if (resp.status_code != 200) {
        throw std::runtime_error("Storage server error when listing items: " + get_response_error_message(resp));
    }

    json j = json::parse(resp.text);
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

StorageItemList StorageClient::list_items(StorageItemTags const& tags, size_t limit) {
    auto query_parameters = tags_to_query_parameters(tags);
    query_parameters.Add({"_limit", std::to_string(limit)});

    auto resp = cpr::Get(cpr::Url(base_url + "/v1/blobs"), query_parameters);
    return read_items_from_response(resp);
}

StorageItemList StorageClient::get_next_page_of_items(StorageItemList const& page) {
    auto resp = cpr::Get(cpr::Url(page.continuation));
    return read_items_from_response(resp);
}

std::shared_ptr<std::istream> StorageClient::get_latest_item(StorageItemTags const& tags) {
    auto resp = cpr::Get(cpr::Url(base_url + "/v1/blobs/data/latest"), tags_to_query_parameters(tags));
    if (resp.status_code == 404) {
        return {};
    }

    if (resp.status_code != 200) {
        throw std::runtime_error("Storage server error when getting latest items: " + get_response_error_message(resp));
    }
    auto stream = std::make_shared<std::stringstream>();
    stream->str(resp.text);
    return stream;
}

std::shared_ptr<std::istream> StorageClient::get_item_by_url(const std::string& url) {
    auto resp = cpr::Get(cpr::Url(url));
    if (resp.status_code == 404) {
        return {};
    }

    if (resp.status_code != 200) {
        throw std::runtime_error("Storage server error when getting item: " + get_response_error_message(resp));
    }

    auto stream = std::make_shared<std::stringstream>();
    stream->str(resp.text);
    return stream;
}

StorageItem StorageClient::store_item(StorageItemTags const& tags, std::istream& data,
                                      std::optional<std::chrono::seconds> time_to_live) {
    auto query_parameters = tags_to_query_parameters(tags);
    if (time_to_live) {
        query_parameters.Add({"_ttl", std::to_string(time_to_live->count()) + "s"});
    }

    std::string s(std::istreambuf_iterator<char>(data), {});
    cpr::Body body(s);
    auto resp = cpr::Post(cpr::Url(base_url + "/v1/blobs/data"), query_parameters, body);
    if (resp.status_code != 201) {
        throw std::runtime_error("Storage server error when storing item: " + get_response_error_message(resp));
    }

    return storage_item_from_json(json::parse(resp.text));
}

std::optional<std::string> StorageClient::health_check() {
    auto response = cpr::Get(cpr::Url(base_url + "/healthcheck"));
    if (response.status_code == 200) {
        return std::nullopt;
    }

    return "Storage server error: " + get_response_error_message(response);
}
} // namespace Gadgetron::Storage