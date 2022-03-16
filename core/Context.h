#pragma once
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include "Storage.h"

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
                storage
            },
            args{args},
            storage_address{storage_address} {}

        Args args;
        StorageAddress storage_address;
    };
}
