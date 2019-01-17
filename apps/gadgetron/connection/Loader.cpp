#include "Loader.h"

#include <memory>

#include "stream/Stream.h"

namespace {
    using namespace Gadgetron::Core;

    using reader_factory = std::unique_ptr<Reader>();
    using writer_factory = std::unique_ptr<Writer>();
}

namespace Gadgetron::Server::Connection {

    Loader::Loader(const Context &context) : context(context) {}

    boost::filesystem::path Loader::make_library_path(const std::string &shared_library_name) const {
        return context.paths.gadgetron_home / "lib" / shared_library_name;
    }

    boost::dll::shared_library Loader::load_library(const std::string &shared_library_name) {

        auto lib = boost::dll::shared_library(
                make_library_path(shared_library_name),
                boost::dll::load_mode::append_decorations
        );

        libraries.push_back(lib);
        return lib;
    }

    std::unique_ptr<Reader> Loader::load(const Config::Reader &conf) {
        auto factory = load_factory<reader_factory>("reader_factory_export_", conf.classname, conf.dll);
        return factory();
    }

    std::unique_ptr<Writer> Loader::load(const Config::Writer &conf) {
        auto factory = load_factory<writer_factory>("writer_factory_export_", conf.classname, conf.dll);
        return factory();
    }

    std::unique_ptr<Stream::Stream> Loader::load(const Config::Stream &conf) {
        return std::make_unique<Stream>(conf, context, *this);
    }
}
