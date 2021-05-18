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
#include <range/v3/view/transform.hpp>
#include <range/v3/range/conversion.hpp>
#include <regex>

namespace beast = boost::beast;     // from <boost/beast.hpp>
namespace http = beast::http;       // from <boost/beast/http.hpp>
namespace net = boost::asio;        // from <boost/asio.hpp>
using tcp = net::ip::tcp;
using json = nlohmann::json;

namespace {
    using namespace Gadgetron::Core;

    constexpr int http_version = 11;

    auto make_json_request(http::verb method, const std::string &host, const std::string &target, const json &j) {

        http::request<http::string_body> req{method, target, http_version};
        req.set(http::field::host, host);
        req.set(http::field::user_agent, BOOST_BEAST_VERSION_STRING);
        req.set(http::field::content_type, "application/json");
        req.body() = j.dump();
        req.prepare_payload();
        return req;
    }

    auto make_empty_request(http::verb method, const std::string &host, const std::string &target) {

        http::request<http::empty_body> req{method, target, http_version};
        req.set(http::field::host, host);
        req.set(http::field::user_agent, BOOST_BEAST_VERSION_STRING);
        return req;
    }

    auto connect_stream(net::io_context &ioc, const std::string &host, const std::string &port) {
        tcp::resolver resolver(ioc);
        beast::tcp_stream stream(ioc);
        auto const results = resolver.resolve(host, port);
        stream.connect(results);
        return stream;
    }

    json get_content(const std::string &server_address, const std::string &port, const std::string &group,
                     const std::string &subject, const std::string &key) {

        net::io_context ioc;
        auto stream = connect_stream(ioc, server_address, port);
        auto json_body = json{{"storagespace", group},
                              {"subject",      subject},
                              {"key",          key}};
        auto req = make_json_request(http::verb::get, server_address, "/v1/data", json_body);
        http::write(stream, req);

        beast::flat_buffer buffer;
        http::response<http::string_body> response;
        http::read(stream, buffer, response);
        stream.socket().shutdown(tcp::socket::shutdown_both);
        return json::parse(response.body());
    }

    std::vector<char> fetch_data(const std::string &server_address, const std::string port, const std::string &path) {
        net::io_context ioc;
        auto stream = connect_stream(ioc, server_address, port);
        auto req = make_empty_request(http::verb::get, server_address, path);
        http::write(stream, req);
        beast::flat_buffer buf;

        http::response_parser<http::vector_body<char>> response_parser;
        response_parser.body_limit(128ull * 1024ull * 1024ull * 1024ull); //We support files up to 128GB. For now.
        http::read(stream, buf, response_parser);
        stream.socket().shutdown(tcp::socket::shutdown_both);
        return std::move(response_parser.get().body());
    }

    json store_request(const std::string &server_address, const std::string &port, const std::string &group,
                       const std::string &subject, const std::string &key,
                       const boost::posix_time::time_duration &duration) {

        net::io_context ioc;
        auto stream = connect_stream(ioc, server_address, port);
        auto json_body = json{{"storagespace", group},
                              {"subject",      subject},
                              {"key",          key},
                              {"storage_duration",     boost::posix_time::to_simple_string(duration)}};
        auto req = make_json_request(http::verb::post, server_address, "/v1/data", json_body);
        http::write(stream, req);

        beast::flat_buffer buffer;
        http::response<http::string_body> response;
        http::read(stream, buffer, response);
        stream.socket().shutdown(tcp::socket::shutdown_both);
        if (to_status_class(response.result()) != http::status_class::successful) throw std::runtime_error("Storage server reported error " + response.body());
        return json::parse(response.body());
    }

    void store_content(const std::string &server_address, const std::string &port, const std::string &path,
                       std::vector<char> data) {
        net::io_context ioc;
        auto stream = connect_stream(ioc, server_address, port);

        http::request<http::vector_body<char>> req{http::verb::patch, path, http_version};
        req.set(http::field::host, server_address);
        req.set(http::field::user_agent, BOOST_BEAST_VERSION_STRING);
        req.body() = std::move(data);
        req.prepare_payload();
        http::write(stream, std::move(req));
        http::response<http::string_body> response;
        beast::flat_buffer buffer;
        http::read(stream, buffer, response);
        stream.socket().shutdown(tcp::socket::shutdown_both);
        if (to_status_class(response.result()) != http::status_class::successful) throw std::runtime_error("Storage server reported error " + response.body());

    }

}

Gadgetron::Storage::RESTStorageClient::RESTStorageClient(const Address& address,
                                                         const std::string &group)
        : server_address(address.host), port(address.port), group(group)  {}


std::vector<std::string> Gadgetron::Storage::RESTStorageClient::content(const std::string& subject, const std::string &key) const {
    auto response = get_content(server_address, port, group, subject, key);

    return response | ranges::views::transform(
            [](const auto &json_value) { return json_value["storage_path"].template get<std::string>(); }) |
           ranges::to<std::vector>();
}

std::vector<char> Gadgetron::Storage::RESTStorageClient::fetch(const std::string &uuid) const {

    return fetch_data(server_address, port, uuid);
}

void Gadgetron::Storage::RESTStorageClient::store(const std::string& subject, const std::string &key, const std::vector<char> &value,
                                                  boost::posix_time::time_duration duration) {
    auto response = store_request(server_address, port, group, subject, key, duration);
    auto blob_path = response["storage_path"];
    store_content(server_address, port, blob_path, value);

}
namespace {
    using namespace Gadgetron::Core;

    struct IDs {
        std::string patientID;
        std::string scannerID;
        std::string studyID;
    };

    optional<IDs> extract_IDS_from_measurement(const std::string &measurementID) {
        std::regex reg("(.*?)_(.*?)_(.*?)_(.*)");
        std::smatch match;

        if (std::regex_match(measurementID, match, reg)) {
            return IDs{match[2], match[1],match[3]};
        }
        return {};
    }


    optional<std::string> scannerID(const ISMRMRD::IsmrmrdHeader &header) {

        if (auto system_info = header.acquisitionSystemInformation)
            if (auto s = system_info->stationName)
                return s.value();

        if (auto measInfo = header.measurementInformation) {
            if (auto id = measInfo->measurementID) {
                if (auto extracted = extract_IDS_from_measurement(id.value())) {
                    return extracted->scannerID;
                }
            }
        }
        return {};
    }


    optional<std::string> patientID(const ISMRMRD::IsmrmrdHeader &header) {
        if (auto subject = header.subjectInformation) {
            if (subject->patientID) {
                return subject->patientID.value();
            }
        }

        if (auto measInfo = header.measurementInformation) {
            if (auto id = measInfo->measurementID) {
                if (auto extracted = extract_IDS_from_measurement(id.value())) {
                    return extracted->patientID;
                }
            }
        }

        return {};

    }
    optional<std::string> studyID(const ISMRMRD::IsmrmrdHeader &header) {
        if (auto subject = header.studyInformation) {
            if (subject->studyID) {
                return subject->studyID.value();
            }
        }

        if (auto measInfo = header.measurementInformation) {
            if (auto id = measInfo->measurementID) {
                if (auto extracted = extract_IDS_from_measurement(id.value())) {
                    return extracted->studyID;
                }
            }
        }

        return {};

    }
    optional<std::string> patientStudyID(const ISMRMRD::IsmrmrdHeader& header){
        auto patient = patientID(header);
        auto study = studyID(header);

        if (patient && study)
            return *patient + "/"  + *study;
        if (study)
            return *study;
        return {};
    }
    optional<std::string> debugID(const ISMRMRD::IsmrmrdHeader &header) {
        return "";
    }


    optional<std::string> measurementID(const ISMRMRD::IsmrmrdHeader& header){
        if (header.measurementInformation)
            if (header.measurementInformation->measurementID)
                return *header.measurementInformation->measurementID;
        return {};
    }


}
Gadgetron::StorageSpaces Gadgetron::Storage::setup_storage(const Address& address,
                                          const ISMRMRD::IsmrmrdHeader &header) {

    auto create_storage = [&](const auto& group, const auto& id, const auto& duration){
        return StorageSpace(std::make_shared<RESTStorageClient>(address, group), id,duration);
    };

        return {
                create_storage("session", patientStudyID(header),boost::posix_time::time_duration(48,0,0)),
                create_storage("scanner", scannerID(header),boost::posix_time::time_duration(48,0,0)),
                MeasurementSpace(std::make_shared<RESTStorageClient>(address, "measurement" ), measurementID(header), boost::posix_time::time_duration(48,0,0))
        };
}
