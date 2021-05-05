#pragma once
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include "Storage.h"

namespace Gadgetron::Core {

    struct Context {


        using Header =
                ISMRMRD::IsmrmrdHeader;

        struct Paths {
            boost::filesystem::path gadgetron_home;
            boost::filesystem::path working_folder;
        };

        Header header;
        Paths  paths;
        Storage storage;
    };

    struct StreamContext : Context {
        using Args = boost::program_options::variables_map;
        StreamContext(ISMRMRD::IsmrmrdHeader header, const Paths paths, Storage storage, const Args args) : Context{std::move(header),paths, std::move(storage)},args{args} {}
        Args   args;
    };


}


