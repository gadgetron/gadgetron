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
            std::vector<std::function<Reader*(void)>> reader_factories;
            std::vector<std::function<Writer*(void)>> writer_factories;
            std::function<std::shared_ptr<Node>(const Header&)> stream_factory;
        };

        struct Header {

        };
    };



    //Should maybe not be here? No worries, cleanups are always fun







}



#endif //GADGETRON_CONTEXT_H
