
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>

#include "log.h"
#include "paths.h"

#include "system_info.h"

#include "Server.h"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace Gadgetron::Server;

int main(int argc, char *argv[]) {

    options_description desc("Allowed options:");
    desc.add_options()
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

    GINFO("Running on port %d\n", args["port"].as<unsigned short>());

    try {
        // Ensure working directory exists.
        create_directories(args["dir"].as<path>());

        boost::asio::io_service io_service;
        Server server(io_service, args);
        io_service.run();
    }
    catch (std::exception &e) {
        GERROR_STREAM(e.what() << std::endl);
        return 1;
    }

    return 0;
}
