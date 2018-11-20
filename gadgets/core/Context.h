#ifndef GADGETRON_CONTEXT_H
#define GADGETRON_CONTEXT_H

#include <boost/filesystem.hpp>

namespace Gadgetron::Core {

    struct Context {

        struct Paths {
            Paths(const boost::filesystem::path &home, const boost::filesystem::path &work);

            const boost::filesystem::path gadgetron_home;
            const boost::filesystem::path working_folder;
        };

        struct Config {

        };

        struct Header {

        };
    };
}



#endif //GADGETRON_CONTEXT_H
