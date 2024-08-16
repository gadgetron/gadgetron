#include "Loader.h"

#include <memory>

#include "nodes/Stream.h"

namespace {
    using namespace Gadgetron::Core;

    using reader_factory = std::unique_ptr<Reader>();
    using writer_factory = std::unique_ptr<Writer>();
}

namespace Gadgetron::Server::Connection {

    Loader::Loader(const StreamContext &context) : context(context) {}

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

    std::unique_ptr<Reader> Loader::load(const Config::Reader &conf) {
        auto factory = load_factory<reader_factory>("reader_factory_export_", conf.classname, conf.dll);
        return factory();
    }

    std::unique_ptr<Writer> Loader::load(const Config::Writer &conf) {
        auto factory = load_factory<writer_factory>("writer_factory_export_", conf.classname, conf.dll);
        return factory();
    }

    std::unique_ptr<Nodes::Stream> Loader::load(const Config::Stream &conf) {
        return std::make_unique<Nodes::Stream>(conf, context, *this);
    }

    std::map<uint16_t, std::unique_ptr<Reader>> Loader::load_readers(const std::vector<Config::Reader> &configs) {

        std::map<uint16_t, std::unique_ptr<Reader>> readers{};

        for (auto &config : configs) {
            auto reader = load(config);
            uint16_t slot = config.slot.value_or(reader->slot());
            readers[slot] = std::move(reader);
        }

        return readers;
    }

    std::vector<std::unique_ptr<Writer>> Loader::load_writers(const std::vector<Config::Writer> &configs) {

        std::vector<std::unique_ptr<Writer>> writers{};

        for (auto &writer_config : configs) {
            writers.emplace_back(load(writer_config));
        }

        return std::move(writers);
    }
}
