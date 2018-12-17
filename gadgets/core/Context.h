#ifndef GADGETRON_CONTEXT_H
#define GADGETRON_CONTEXT_H

#include <boost/filesystem.hpp>
#include <ismrmrd/ismrmrd.h>

#include <ismrmrd/xml.h>

namespace Gadgetron::Core {

    struct Context {

        using Header =
                ISMRMRD::IsmrmrdHeader;

        struct Paths {
            Paths(const boost::filesystem::path &home, const boost::filesystem::path &work);

            const boost::filesystem::path gadgetron_home;
            const boost::filesystem::path working_folder;
        };

        Header header;
        Paths  paths;
    };
}


#endif //GADGETRON_CONTEXT_H
