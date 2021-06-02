
#include <boost/program_options.hpp>
#include <tuple>
#include <string>
#include <optional>

#include "StorageServer.h"
#include "gadgetron_paths.h"
#include "storage.h"
#include "log.h"

using namespace Gadgetron::Storage;

using namespace boost::filesystem;
using namespace boost::program_options;


namespace Gadgetron::Server {

    std::tuple<std::string, std::optional<StorageServer>> ensure_storage_server(const variables_map& args) {

        if (args.count("storage_address")) {
            return {std::string(args["storage_address"].as<std::string>()), std::nullopt};
        }

        auto server = StorageServer(
            args["storage_port"].as<unsigned short>(),
            args["database_dir"].as<path>(),
            args["storage_dir"].as<path>()
        );

        GINFO_STREAM("Running storage server on port " + std::to_string(server.port()))

        return {"http://localhost:" + std::to_string(server.port()), std::move(server)};
    }
}