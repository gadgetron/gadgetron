#ifndef GADGETRON_CONTEXT_H
#define GADGETRON_CONTEXT_H

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

        Header  header;
        Paths   paths;
        Storage storage;
    };

}


#endif //GADGETRON_CONTEXT_H
