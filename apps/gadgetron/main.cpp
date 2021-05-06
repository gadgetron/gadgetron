
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>

#include "log.h"
#include "paths.h"
#include "initialization.h"

#include "system_info.h"
#include "gadgetron_config.h"

#include "Server.h"
#include "StorageServer.h"

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
            ("storage_port,s", value<unsigned short>()->default_value(0), "Port on which to run the Storage server")
            ("database_dir,D",value<path>()->default_value(default_database_folder()),"Directory in which to store the database")
            ("storage_dir,S", value<path>()->default_value(default_storage_folder()),"Directory in which to store data blob");
    options_description desc;
    desc.add(gadgetron_options).add(storage_options);
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

    try {
        check_environment_variables();
        configure_blas_libraries();
        set_locale();

        // Ensure working directory exists.
        create_directories(args["dir"].as<path>());

        auto sessionsServer = Gadgetron::Storage::StorageServer(args["storage_port"].as<unsigned short>(), args["database_dir"].as<path>(), args["storage_dir"].as<path>());

        Server server(args, {"localhost", std::to_string(sessionsServer.port())});
        server.serve();
    }
    catch (std::exception &e) {
        GERROR_STREAM(e.what() << std::endl);
        return 1;
    }

    return 0;
}
