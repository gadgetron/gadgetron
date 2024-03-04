
#include <optional>
#include <regex>
#include <string>
#include <tuple>
#include <utility>
#include <thread>

#include <nlohmann/json.hpp>

#include "IsmrmrdContextVariables.h"
#include "Process.h"
#include "gadgetron_paths.h"
#include "log.h"

#include "storage.h"

using namespace Gadgetron::Storage;
using namespace boost::filesystem;
using namespace boost::program_options;

namespace Gadgetron::Server {

void invoke_storage_server_health_check(std::string base_address) {
    GINFO_STREAM("Verifying connectivity to storage server...")
    StorageClient c(base_address);
    for (int i = 0;; i++) {
        auto err = c.health_check();
        if (!err.has_value()) {
            GINFO_STREAM("Received successful response from storage server.");
            break;
        }

        if (i == 49) {
            throw std::runtime_error(*err);
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
}

std::tuple<std::string, std::optional<boost::process::child>> ensure_storage_server(const variables_map& args) {
    if (args.count("disable_storage") && args["disable_storage"].as<bool>()) {
        return {"", std::nullopt};
    }

    if (args.count("storage_address")) {
        auto uri = args["storage_address"].as<std::string>();
        invoke_storage_server_health_check(uri);
        return {uri, std::nullopt};
    }

    auto port = args["storage_port"].empty() ? args["port"].as<unsigned short>() + 110
                                             : args["storage_port"].as<unsigned short>();
    auto environment = boost::this_process::environment();
    environment.set("MRD_STORAGE_SERVER_PORT", std::to_string(port));
    environment.set("MRD_STORAGE_SERVER_DATABASE_CONNECTION_STRING",
                    (args["database_dir"].as<path>() / "metadata.db").string());
    environment.set("MRD_STORAGE_SERVER_STORAGE_CONNECTION_STRING", (args["storage_dir"].as<path>()).string());

    auto storage_executable = boost::process::search_path("mrd-storage-server");
    if (storage_executable.empty()) {
        throw std::runtime_error("Failed to find MRD Storage Server.\n"
                                 "Please ensure 'mrd-storage-server' is found on your PATH.\n"
                                 "See: https://github.com/ismrmrd/mrd-storage-server/");
    }

    GDEBUG_STREAM("Found storage server: " + storage_executable.string())
    GINFO_STREAM("Starting storage server on port " + environment.get("MRD_STORAGE_SERVER_PORT"))

    auto uri = "http://localhost:" + environment.get("MRD_STORAGE_SERVER_PORT");
    auto process =
        Process::child(storage_executable, "--require-parent-pid", // have child process exit if parent crashes
                       std::to_string(boost::this_process::get_id()), 
                       boost::process::std_out > boost::process::null, boost::process::std_err > stderr, environment);

    invoke_storage_server_health_check(uri);

    return {uri, std::move(process)};
}

StorageSpaces setup_storage_spaces(const std::string& address, const ISMRMRD::IsmrmrdHeader& header) {
    auto client = std::make_shared<StorageClient>(address);
    IsmrmrdContextVariables variables(header);
    auto ttl = std::chrono::hours(48);

    return {std::make_shared<SessionSpace>(client, variables, ttl),
            std::make_shared<ScannerSpace>(client, variables, ttl),
            std::make_shared<MeasurementSpace>(client, variables, ttl)};
}
} // namespace Gadgetron::Server
