//
// Created by dch on 5/21/21.
//
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "gadgetron_paths.h"
#include <iostream>
#include "StorageServer.h"

using namespace boost::program_options;
using namespace boost::filesystem;
using namespace Gadgetron;

int main(int argc, char** argv){


    options_description storage_options("Storage options");
    storage_options.add_options()
            ("storage_port,s", value<unsigned short>()->default_value(9112), "Port on which to run the Storage server")
            ("database_dir,D",value<path>()->default_value(Server::default_database_folder()),"Directory in which to store the database")
            ("storage_dir,S", value<path>()->default_value(Server::default_storage_folder()),"Directory in which to store data blob");
    options_description desc;
    desc.add(storage_options);
    variables_map args;
    store(parse_command_line(argc, argv, desc), args);
    notify(args);
    if (args.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    auto sessionsServer = Gadgetron::Storage::StorageServer(args["storage_port"].as<unsigned short>(), args["database_dir"].as<path>(), args["storage_dir"].as<path>());
    std::cout << "Running Gadgetron Storage server on port " << sessionsServer.port() << std::endl;
    sessionsServer.run_forever();
}