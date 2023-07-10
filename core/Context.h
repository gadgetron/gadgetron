#pragma once
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include "StorageSetup.h"

inline std::map<std::string, std::string> GetParameters(const boost::program_options::variables_map& args) {
    std::map<std::string, std::string> parameters;
    if (args.count("parameter")) {
        auto params = args["parameter"].as<std::vector<std::pair<std::string, std::string>>>();
        for (auto &arg : params) {
            parameters[arg.first] = arg.second;
        }
    } 
    return parameters;
}

namespace Gadgetron::Core {

    struct Context {
        using Header = ISMRMRD::IsmrmrdHeader;

        struct Paths {
            boost::filesystem::path gadgetron_home;
            boost::filesystem::path working_folder;
        };

        Header header;
        Paths  paths;
        StorageSpaces storage;
        std::map<std::string, std::string> parameters;
    };

    struct StreamContext : Context {
        using Args = boost::program_options::variables_map;
        using StorageAddress = std::string;

        StreamContext(
            ISMRMRD::IsmrmrdHeader header,
            const Paths paths,
            const Args args,
            const StorageAddress storage_address,
            StorageSpaces storage
        ) : Context{
                std::move(header),
                paths,
                storage,
                GetParameters(args)
            },
            args{args},
            storage_address{storage_address} {}

        Args args;
        StorageAddress storage_address;
    };
}
