
#include <regex>
#include <tuple>
#include <string>
#include <utility>
#include <optional>

#include <boost/process.hpp>
#include <boost/program_options.hpp>

#include <cpr/cpr.h>
#include <nlohmann/json.hpp>

#include "Process.h"
#include "gadgetron_paths.h"
#include "log.h"

#include "storage.h"


namespace {
    using json = nlohmann::json;

    cpr::Parameters as_parameters(const std::map<std::string, std::string> &query) {
        auto parameters = cpr::Parameters{};
        for (auto const& [key, value] : query) parameters.Add({key, value});
        return parameters;
    }

    std::vector<char> as_data(const std::string &str) {
        std::vector<char> data(str.length());
        std::copy(str.begin(), str.end(), data.begin());
        return data;
    }

    class MRDStorageClient {
      private:
        const std::string address;

      public:
        explicit MRDStorageClient(std::string address) : address{std::move(address)} {}

        json get(const std::string &uri) {
            auto response = cpr::Get(cpr::Url(this->address));
            if (response.status_code != 200) throw std::runtime_error(response.status_line);
            return json::parse(response.text);
        }

        json search(const std::map<std::string, std::string> &query) {
            auto response = cpr::Get(cpr::Url(this->address + "/v1/blobs"), as_parameters(query));
            if (response.status_code != 200) throw std::runtime_error(response.status_line);
            return json::parse(response.text);
        }

        std::vector<char> fetch(const std::string &uri) {
            auto response = cpr::Get(cpr::Url(uri));
            if (response.status_code != 200) throw std::runtime_error(response.status_line);
            return as_data(response.text);
        }

        json store(const std::vector<char> &data,
                   const std::map<std::string, std::string> &query) {
            auto response = cpr::Post(cpr::Url(this->address + "/v1/blobs/data"), as_parameters(query),
                                      cpr::Body(data.data(), data.size()));
            if (response.status_code != 201) throw std::runtime_error(response.status_line);
            return json::parse(response.text);
        }
    };
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

        if (patient && study) {
            std::string result =  *patient + "/" + *study;
            return result;

        }

        if (study) return *study;
        return {};
    }

    optional<std::string> measurementID(const ISMRMRD::IsmrmrdHeader& header){
        if (header.measurementInformation)
            if (header.measurementInformation->measurementID)
                return *header.measurementInformation->measurementID;
        return {};
    }
}

using namespace Gadgetron::Storage;
using namespace boost::filesystem;
using namespace boost::program_options;

namespace Gadgetron::Server {

    std::tuple<std::string, std::optional<boost::process::child>>
    ensure_storage_server(const variables_map& args) {

        if (args.count("storage_address")) {
            return {args["storage_address"].as<std::string>(), std::nullopt};
        }

        auto environment = boost::this_process::environment();
        environment.set("MRD_STORAGE_SERVER_PORT", std::to_string(args["storage_port"].as<uint16_t>()));
        environment.set("MRD_STORAGE_SERVER_DATABASE_CONNECTION_STRING",
                        (args["database_dir"].as<path>() / "metadata.db").string());
        environment.set("MRD_STORAGE_SERVER_STORAGE_CONNECTION_STRING",
                        (args["storage_dir"].as<path>()).string());

        auto storage_executable = boost::process::search_path("mrd-storage-server");
        if (storage_executable.empty()) {
            throw std::runtime_error("Failed to find MRD Storage Server.\n"
                                     "Please ensure 'mrd-storage-server' is found on your PATH.\n"
                                     "See: https://github.com/ismrmrd/mrd-storage-server/");
        }

        GDEBUG_STREAM("Found storage server: " + storage_executable.string())
        GINFO_STREAM("Starting storage server on port " + environment.get("MRD_STORAGE_SERVER_PORT"))
        return {
            "http://localhost:" + environment.get("MRD_STORAGE_SERVER_PORT"),
            Process::child(storage_executable,
                           boost::process::std_out > boost::process::null,
                           boost::process::std_err > stderr,
                           environment)
        };
    }

    StorageSpaces setup_storage_spaces(
        const std::string& address,
        const ISMRMRD::IsmrmrdHeader &header
    ) {
        GDEBUG_STREAM("Using storage address: " << address);
        return {
            StorageSpace(std::make_shared<StorageClient>(address, "session"), patientStudyID(header), boost::posix_time::time_duration(48,0,0)),
            StorageSpace(std::make_shared<StorageClient>(address, "scanner"), scannerID(header), boost::posix_time::time_duration(48,0,0)),
            MeasurementSpace(std::make_shared<StorageClient>(address, "measurement"), measurementID(header), boost::posix_time::time_duration(48,0,0))
        };
    }
}

namespace {
    std::string format_duration(boost::posix_time::time_duration duration) {
        // TODO: Call std::format (or something) in stead, once C++20 is widely supported.
        return std::to_string(duration.total_seconds()) + "s";
    }
}

namespace Gadgetron::Server {

    StorageClient::StorageClient(
        std::string address,
        std::string group
    ) : address{std::move(address)}, space{std::move(group)} {}

    std::vector<std::string> StorageClient::content(
        const std::string& identifier,
        const std::string& name
    ) const {
        std::vector<std::string> uris;

        auto client = MRDStorageClient(this->address);
        auto current_results = client.search({{"name", name},
                                              {"subject", "$null"},
                                              {this->space, identifier}});

        auto extract_results = [&](const auto &results, const auto &recurse) -> void {
            for (const auto &item : results["items"]) uris.push_back(item["data"]);
            if (results.contains("next")) recurse(client.get(results["next"]), recurse);
        };

        extract_results(current_results, extract_results);
        return uris;
    }

    std::vector<char> StorageClient::fetch(
        const std::string& uri
    ) const {
        auto client = MRDStorageClient(this->address);
        return client.fetch(uri);
    }

    void StorageClient::store(
        const std::string& identifier,
        const std::string& name,
        const std::vector<char>& value,
        boost::posix_time::time_duration duration
    ) {
        auto client = MRDStorageClient(this->address);
        client.store(value,
                     {{"name", name},
                      {"subject", "$null"},
                      {"ttl", format_duration(duration)},
                      {this->space, identifier}});
    }
}