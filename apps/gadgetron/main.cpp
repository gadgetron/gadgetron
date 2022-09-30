#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "log.h"
#include "gadgetron_paths.h"
#include "initialization.h"
#include "storage.h"

#include "system_info.h"
#include "gadgetron_config.h"

#include "Server.h"
#include "StreamConsumer.h"


using namespace boost::filesystem;
using namespace boost::program_options;
using namespace Gadgetron::Server;

int main(int argc, char *argv[]) {
    options_description gadgetron_options("Allowed options:");
    gadgetron_options.add_options()
            ("help,h", "Prints this help message.")
            ("info", "Prints build info about the Gadgetron.")
            ("dir,W",
                value<path>()->default_value(default_working_folder()),
                "Set the Gadgetron working directory.")
            ("home,G",
                value<path>()->default_value(default_gadgetron_home()),
                "Set the Gadgetron home directory.")
            ("port,p",
                value<unsigned short>()->default_value(9002),
                "Listen for incoming connections on this port.")
            ("from_stream, s",
                "Perform reconstruction from a local data stream")
            ("input_path,i",
                value<std::string>(),
                "Input file for binary data to perform a local reconstruction with")
            ("output_path,o",
                value<std::string>(),
                "Output file for binary data as a result of a local reconstruction")
            ("config_name,c",
                value<std::string>(),
                "Filename of the desired gadgetron reconstruction config.");

    options_description storage_options("Storage options");
    storage_options.add_options()
            ("storage_address,E",
                value<std::string>(),
                "External address of a storage server. If not provided, a storage server will be started.")
            ("storage_port,s",
                value<unsigned short>(),
                "Port on which to run the storage server. "
                "If no port is provided, a port offset from the port argument (-p) is selected.")
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
            GINFO_STREAM(desc);
            return 0;
        }

        if (args.count("info")) {
            std::stringstream str;
            Info::print_system_information(str);
            GINFO(str.str().c_str());
            return 0;
        }

        GINFO("Gadgetron %s [%s]\n", GADGETRON_VERSION_STRING, GADGETRON_GIT_SHA1_HASH);

        // Ensure working directory exists.
        create_directories(args["dir"].as<path>());

        auto [storage_address, storage_server] = ensure_storage_server(args);

        if(!args.count("from_stream"))
        {
            GINFO("Running on port %d\n", args["port"].as<unsigned short>());
            Server server(args, storage_address);
            server.serve();
        }
        else
        {
            auto cfg = args["config_name"].as<std::string>();
            StreamConsumer consumer(args, storage_address);

            if(args.count("input_path") && args.count("output_path"))
            {
                auto input_stream = std::ifstream(args["input_path"].as<std::string>());
                auto output_stream = std::ofstream(args["output_path"].as<std::string>());
                consumer.consume(input_stream, output_stream, cfg);
                output_stream.close();
            }
            else if(args.count("input_path"))
            {
                auto input_stream = std::ifstream(args["input_path"].as<std::string>());
                consumer.consume(input_stream, std::cout, cfg);
                std::flush(std::cout);
            }
            else if(args.count("output_path"))
            {
                auto output_stream = std::ofstream(args["output_path"].as<std::string>());
                consumer.consume(std::cin, output_stream, cfg);
                output_stream.close();
            }
            else
            {
                consumer.consume(std::cin, std::cout, cfg);
                std::flush(std::cout);
            }
        }
    }
    catch (std::exception &e)
    {
        GERROR_STREAM(e.what() << std::endl);
        std::exit(EXIT_FAILURE);
    }
    catch(...)
    {
        GERROR_STREAM("Unhandled exception, exiting" << std::endl);
        std::exit(EXIT_FAILURE);
    }

    return 0;
}
