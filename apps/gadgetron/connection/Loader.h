#pragma once

#include <map>
#include <memory>

#include "stream/Stream.h"
#include "Config.h"

#include "Context.h"
#include "Reader.h"
#include "Writer.h"

namespace Gadgetron::Server::Connection::Stream {
    class Stream;
}

namespace Gadgetron::Server::Connection {

    class Loader {
        using Context = Core::Context;
        using Reader  = Gadgetron::Core::Reader;
        using Writer  = Gadgetron::Core::Writer;

        using GadgetProperties = Core::GadgetProperties;
    public:
        explicit Loader(const Context &);

        std::unique_ptr<Reader> load(const Config::Reader &);
        std::unique_ptr<Writer> load(const Config::Writer &);
        std::unique_ptr<Stream::Stream> load(const Config::Stream &);

        template<class RESULT>
        using generic_factory = std::unique_ptr<RESULT>(
                const Context &,
                const GadgetProperties &
        );

        template<class FACTORY>
        FACTORY& load_factory(const std::string &prefix, const std::string &classname, const std::string &dll) {
            auto library = load_library(dll);
            return library.get_alias<FACTORY>(prefix + classname);
        }

        template<class CONFIG>
        std::map<uint16_t, std::unique_ptr<Reader>> load_readers(CONFIG config) {

            std::map<uint16_t, std::unique_ptr<Reader>> readers{};

            for (auto &reader_config : config.readers) {
                auto reader = load(reader_config);
                uint16_t slot = reader_config.slot.value_or(reader->slot());
                readers[slot] = std::move(reader);
            }

            return std::move(readers);
        }

        template<class CONFIG>
        std::vector<std::unique_ptr<Writer>> load_writers(CONFIG config) {

            std::vector<std::unique_ptr<Writer>> writers{};

            for (auto &writer_config : config.writers) {
                writers.emplace_back(load(writer_config));
            }

            return std::move(writers);
        }


    private:
        boost::dll::shared_library load_library(const std::string &shared_library_name);
        boost::filesystem::path make_library_path(const std::string &shared_library_name) const;

        const Context context;

        std::vector<boost::dll::shared_library> libraries = std::vector<boost::dll::shared_library>();
    };
}

