
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>

#include "log.h"
#include "paths.h"

#include "gadgetron_config.h"
#include "system_info.h"

#include "Server.h"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace Gadgetron::Server;

int main(int argc, char* argv[]) {

    Settings settings;

    options_description desc("Allowed options:");
    desc.add_options()("help", "Prints this help message.")("info", "Prints build info about the Gadgetron.")("dir,W",
        value<path>(&settings.paths.working_folder)->default_value(default_working_folder()),
        "Set the Gadgetron working directory.")("home,G",
        value<path>(&settings.paths.gadgetron_home)->default_value(default_gadgetron_home()),
        "Set the Gadgetron home directory.")("port,p", value<unsigned short>(&settings.port)->default_value(9002),
        "Listen for incoming connections on this port.")("storage_server,S", value<std::string>(),
        "Address for the Session storage server. Will start a local one if not specified");

    variables_map args;
    store(parse_command_line(argc, argv, desc), args);
    notify(args);

    if (args.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    if (args.count("info")) {
        Info::print_system_information(std::cout);
        return 0;
    }

    GINFO("Gadgetron %s [%s]\n", GADGETRON_VERSION_STRING, GADGETRON_GIT_SHA1_HASH);
    GINFO("Running on port %d\n", args["port"].as<unsigned short>());
    auto storage_process = std::optional<boost::process::child>{};
    try {
        // Ensure working directory exists.
        create_directories(args["dir"].as<path>());

        if (!args.count("storage_server")) {
            auto [proc, address]     = start_storage_server(settings.paths.working_folder);
            storage_process          = std::move(proc);
            settings.storage_address = address;
        } else {
            settings.storage_address = args["storage_server"].as<std::string>();
        }

        Server server(settings);
        server.serve();
    } catch (std::exception& e) {
        GERROR_STREAM(e.what() << std::endl);
        return 1;
    }

    return 0;
}
