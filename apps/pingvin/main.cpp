#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "log.h"
#include "initialization.h"

#include "system_info.h"
#include "pingvin_config.h"

#include "StreamConsumer.h"


using namespace boost::filesystem;
using namespace boost::program_options;
using namespace Gadgetron::Main;

using gadget_parameter = std::pair<std::string, std::string>;

std::istream& operator>>(std::istream& in, gadget_parameter& param) {
    std::string token;
    in >> token;
    // parse <key>=<value> into a gadget_parameter
    auto pos = token.find('=');
    if (pos == std::string::npos) {
        throw std::runtime_error("Invalid gadget parameter: " + token);
    }
    param.first = token.substr(0, pos);
    param.second = token.substr(pos + 1);
    return in;
}

std::ostream& operator<<(std::ostream& out, const gadget_parameter& param) {
    out << param.first << "=" << param.second;
    return out;
}

int main(int argc, char *argv[]) {
    options_description gadgetron_options("Allowed options:");
    gadgetron_options.add_options()
            ("help,h", "Prints this help message.")
            ("info", "Prints build info about Pingvin.")
            ("home,G",
                value<path>()->default_value(Info::default_pingvin_home()),
                "Set the Pingvin home directory.")
            ("input,i",
                value<std::string>(),
                "Input file for binary data to perform a local reconstruction with")
            ("output,o",
                value<std::string>(),
                "Output file for binary data as a result of a local reconstruction")
            ("config,c",
                value<std::string>(),
                "Filename of the desired Pingvin reconstruction config.")
            ("parameter",
                value<std::vector<gadget_parameter>>(),
                "Parameter to be passed to the Pingvin reconstruction config. Multiple parameters can be passed."
                "Format: --parameter <name>=<value> --parameter <name>=<value> ...");

    options_description desc;
    desc.add(gadgetron_options);

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

        GINFO("Pingvin %s [%s]\n", PINGVIN_VERSION_STRING, PINGVIN_GIT_SHA1_HASH);

        if (!args.count("config"))
        {
            GERROR_STREAM("No config file provided. Use --config/-c");
            return 1;
        }

        auto cfg = args["config"].as<std::string>();
        StreamConsumer consumer(args);

        std::unique_ptr<std::istream> input_file;
        if (args.count("input")) {
            input_file = std::make_unique<std::ifstream>(args["input"].as<std::string>());
            if (!input_file->good()) {
                GERROR_STREAM("Could not open input file: " << args["input"].as<std::string>());
                return 1;
            }
        }

        std::unique_ptr<std::ostream> output_file;
        if (args.count("output")) {
            output_file = std::make_unique<std::ofstream>(args["output"].as<std::string>());
            if (!output_file->good()) {
                GERROR_STREAM("Could not open output file: " << args["output"].as<std::string>());
                return 1;
            }
        }

        std::istream& input = input_file ? *input_file : std::cin;
        std::ostream& output = output_file ? *output_file : std::cout;
        consumer.consume(input, output, cfg);
        std::flush(output);

        GDEBUG_STREAM("Finished consuming stream");
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
