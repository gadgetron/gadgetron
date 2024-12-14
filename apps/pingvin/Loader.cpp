#include "Loader.h"

#include <memory>

#include "nodes/Stream.h"

namespace Gadgetron::Main {

    Loader::Loader(const Core::StreamContext &context) : context(context) {}

    boost::dll::shared_library Loader::load_library(const std::string &shared_library_name) {
        try
        {
            auto lib = boost::dll::shared_library(
                    shared_library_name,
                    boost::dll::load_mode::append_decorations |
                    boost::dll::load_mode::rtld_global |
                    boost::dll::load_mode::search_system_folders
            );
            libraries.push_back(lib);
            return lib;
        }
        catch( const std::exception & ex )
        {
            std::cerr << ex.what() << std::endl;
            throw;
        }
    }

    std::unique_ptr<Nodes::Stream> Loader::load(const Config::Stream &conf) {
        return std::make_unique<Nodes::Stream>(conf, context, *this);
    }
}
