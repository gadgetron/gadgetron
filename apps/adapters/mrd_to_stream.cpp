
#include <string>
#include <iostream>

#include <boost/asio.hpp>
#include <boost/program_options.hpp>

#include "ismrmrd_msg_adapter.h"

namespace po = boost::program_options;

namespace
{

struct CommandlineArgs
{
    std::string input_path;
    std::string output_path;
    std::string dataset_group;
};

CommandlineArgs ParseCommandLineArgs(int argc, char** argv)
{
    CommandlineArgs args;
    po::variables_map vm;

    po::options_description desc(
        "\n"
        "Converts an ISMRMRD HDF5 file to a binary file of gadgetron messages\n\n"
        "Usage: " + std::string(argv[0]) + " [options]\n\n"
        "Options:");

    desc.add_options()
        ("help,h", "Produce help message")
        (
            "input-path,i",
            po::value<std::string>(&args.input_path)->value_name("Input path")->required(),
            "Path to ISMRMRD HDF5 file being converted."
        )
        (
            "output-path,o",
            po::value<std::string>(&args.output_path)->value_name("Output path"),
            "Path to output binary message data of the converted HDF5 file. Note: If not provided data is written to stdout"
        )
        (
            "dataset-group,d",
            po::value<std::string>(&args.dataset_group)->value_name("Dataset group")->default_value("/dataset"),
            "Dataset within the ISMRMRD file to convert."
        );

    try
    {
        po::store(po::command_line_parser(
                    argc, argv)
                    .options(desc)
                    .run(),
                vm);
        po::notify(vm);
    }
    catch (po::error_with_option_name& e)
    {
        std::cerr << std::endl << e.what() << std::endl << desc << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (vm.count("help"))
    {
        std::cerr << desc << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    return args;
}
}

int main(int argc, char** argv)
{
    try
    {
        CommandlineArgs args = ParseCommandLineArgs(argc, argv);
        auto adapter = IsmrmrdMsgAdapter();

        if(args.output_path.size() > 0)
        {
            std::ofstream stream(args.output_path);
            adapter.convert(stream, args.input_path, args.dataset_group);
            stream.close();
        }
        else
        {
            adapter.convert(std::cout, args.input_path, args.dataset_group);
            close(STDOUT_FILENO);
        }
    }
    catch(const std::exception& exc)
    {
        std::cerr << "Conversion failed, reason: " << exc.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    catch(...)
    {
        std::cerr << "Conversion failed, unhandled exception..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return 0;
}
