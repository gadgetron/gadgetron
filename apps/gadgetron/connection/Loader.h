#pragma once

#include <map>
#include <memory>

#include "config/Config.h"
#include "nodes/Stream.h"

#include "Context.h"
#include "Reader.h"
#include "Writer.h"

namespace Gadgetron::Server::Connection::Nodes {
    class Stream;
}

namespace Gadgetron::Server::Connection {

    class Loader {
        using Context = Core::Context;
        using Reader  = Gadgetron::Core::Reader;
        using Writer  = Gadgetron::Core::Writer;
        using Stream  = Gadgetron::Server::Connection::Nodes::Stream;

        using GadgetProperties = Core::GadgetProperties;
    public:
        explicit Loader(const Core::StreamContext &);

        std::unique_ptr<Reader> load(const Config::Reader &);
        std::unique_ptr<Writer> load(const Config::Writer &);
        std::unique_ptr<Stream> load(const Config::Stream &);

        template<class RESULT>
        using generic_factory = std::unique_ptr<RESULT>(
                const Context &,
                const GadgetProperties &
        );

        template<class FACTORY>
        FACTORY& load_factory(const std::string &prefix, const std::string &classname, const std::string &dll) {
            GINFO_STREAM("loading " << prefix << " - " << classname << " from the dll " << dll);
            auto library = load_library(dll);
            return library.get_alias<FACTORY>(prefix + classname);
        }

        std::map<uint16_t, std::unique_ptr<Reader>> load_readers(const std::vector<Config::Reader> &);
        std::vector<std::unique_ptr<Writer>> load_writers(const std::vector<Config::Writer> &);

        template<class CONFIG>
        std::map<uint16_t, std::unique_ptr<Reader>> load_readers(CONFIG config) {
            return load_readers(config.readers);
        }

        template<class CONFIG>
        std::map<uint16_t, std::unique_ptr<Reader>> load_default_or_custom_readers(CONFIG config) {

            static const std::vector<Config::Reader> default_readers{
                    Config::Reader { "gadgetron_core_readers", "AcquisitionReader", Core::none },
                    Config::Reader { "gadgetron_core_readers", "WaveformReader", Core::none },
                    Config::Reader { "gadgetron_core_readers", "ImageReader", Core::none },
                    Config::Reader { "gadgetron_core_readers", "BufferReader", Core::none },
                    Config::Reader { "gadgetron_core_readers", "IsmrmrdImageArrayReader", Core::none },
                    Config::Reader { "gadgetron_core_readers", "AcquisitionBucketReader", Core::none }
            };

            if (config.readers.empty())
                return load_readers(default_readers);
            return load_readers(config.readers);
        }

        template<class CONFIG>
        std::vector<std::unique_ptr<Writer>> load_writers(CONFIG config) {
            return load_writers(config.writers);
        }

        template<class CONFIG>
        std::vector<std::unique_ptr<Writer>> load_default_or_custom_writers(CONFIG config) {

            static const std::vector<Config::Writer> default_writers{
                    Config::Writer { "gadgetron_core_writers", "AcquisitionWriter" },
                    Config::Writer { "gadgetron_core_writers", "WaveformWriter" },
                    Config::Writer { "gadgetron_core_writers", "ImageWriter" },
                    Config::Writer { "gadgetron_core_writers", "BufferWriter" },
                    Config::Writer { "gadgetron_core_writers", "IsmrmrdImageArrayWriter" },
                    Config::Writer { "gadgetron_core_writers", "AcquisitionBucketWriter" }
            };

            if (config.writers.empty())
                return load_writers(default_writers);
            return load_writers(config.writers);
        }

    private:
        boost::dll::shared_library load_library(const std::string &shared_library_name);

        const Core::StreamContext context;

        std::vector<boost::dll::shared_library> libraries = std::vector<boost::dll::shared_library>();
    };
}

