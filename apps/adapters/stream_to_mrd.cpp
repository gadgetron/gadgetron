#include <string>
#include <iostream>

#include <boost/asio.hpp>
#include <boost/program_options.hpp>

#include "img_msg_adapter.h"

namespace po = boost::program_options;

namespace
{

struct CommandlineArgs
{
    std::string input_path;
    std::string output_path;
    std::string dataset_group;
};

CommandlineArgs parse_command_line_args(int argc, char** argv)
{
    CommandlineArgs args;
    po::variables_map vm;

    po::options_description desc(
        "\n"
        "Converts a binary file of gadgetron messages into an ISMRMRD HDF5 file\n\n"
        "Usage: " + std::string(argv[0]) + " [options]\n\n"
        "Options:");

    desc.add_options()
        ("help,h", "Produce help message")
        (
            "input-path,i",
            po::value<std::string>(&args.input_path)->value_name("Input path"),
            "Path to gadgetron binary data to convert to HDF5. Note: If not provided data is gathered via stdin."
        )
        (
            "output-path,o",
            po::value<std::string>(&args.output_path)->value_name("Output path")->required(),
            "Path to output HDF5 formatted data."
        )
        (
            "dataset-group,g",
            po::value<std::string>(&args.dataset_group)->value_name("Dataset group")->default_value("/dataset"),
            "Dataset to store output images."
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

} // namespace

int main(int argc, char** argv)
{
    try
    {
        CommandlineArgs args = parse_command_line_args(argc, argv);
        auto adapter = ImgMsgAdapter();

        if(args.input_path.size() > 0)
        {
            std::ifstream stream(args.input_path);
            adapter.convert(stream, args.output_path, args.dataset_group);
        }
        else
        {
            adapter.convert(std::cin, args.output_path, args.dataset_group);
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
