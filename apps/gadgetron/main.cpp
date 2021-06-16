
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>

#include "log.h"
#include "gadgetron_paths.h"
#include "initialization.h"
#include "storage.h"

#include "system_info.h"
#include "gadgetron_config.h"

#include "Server.h"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace Gadgetron::Server;

int main(int argc, char *argv[]) {

    options_description gadgetron_options("Allowed options:");
    gadgetron_options.add_options()
            ("help", "Prints this help message.")
            ("info", "Prints build info about the Gadgetron.")
            ("dir,W",
                value<path>()->default_value(default_working_folder()),
                "Set the Gadgetron working directory.")
            ("home,G",
                value<path>()->default_value(default_gadgetron_home()),
                "Set the Gadgetron home directory.")
            ("port,p",
                value<unsigned short>()->default_value(9002),
                "Listen for incoming connections on this port.");

    options_description storage_options("Storage options");
    storage_options.add_options()
            ("storage_address,E",
                value<std::string>(),
                "External address of a storage server. If not provided, a storage server will be started.")
            ("storage_port,s",
                value<unsigned short>()->default_value(9112),
                "Port on which to run the storage server.")
            ("database_dir,D",
                value<path>()->default_value(default_database_folder()),
                "Directory in which to store the storage server database.")
            ("storage_dir,S",
                value<path>()->default_value(default_storage_folder()),
                "Directory in which to store data blobs.");

    options_description desc;
    desc
        .add(gadgetron_options)
        .add(storage_options);

    variables_map args;
    store(parse_command_line(argc, argv, desc), args);
    notify(args);

    try {
        check_environment_variables();
        configure_blas_libraries();
        set_locale();

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

        // Ensure working directory exists.
        create_directories(args["dir"].as<path>());

        auto [storage_address, storage_server] = ensure_storage_server(args);

        Server server(args, storage_address);
        server.serve();
    }
    catch (std::exception &e) {
        GERROR_STREAM(e.what() << std::endl);
        return 1;
    }

    return 0;
}
