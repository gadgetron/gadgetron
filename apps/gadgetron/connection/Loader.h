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
        explicit Loader(const Core::StreamContext &);

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

        std::map<uint16_t, std::unique_ptr<Reader>> load_readers(const std::vector<Config::Reader> &);
        std::vector<std::unique_ptr<Writer>> load_writers(const std::vector<Config::Writer> &);

        template<class CONFIG>
        std::map<uint16_t, std::unique_ptr<Reader>> load_readers(CONFIG config) {
            return load_readers(config.readers);
        }

        template<class CONFIG>
        std::map<uint16_t, std::unique_ptr<Reader>> load_default_and_additional_readers(CONFIG config) {

            static const std::vector<Config::Reader> default_readers{
                    Config::Reader { "gadgetron_core_readers", "AcquisitionReader", boost::none },
                    Config::Reader { "gadgetron_core_readers", "WaveformReader", boost::none },
                    Config::Reader { "gadgetron_core_readers", "ImageReader", boost::none },
                    Config::Reader { "gadgetron_core_readers", "BufferReader", boost::none },
                    Config::Reader { "gadgetron_core_readers", "IsmrmrdImageArrayReader", boost::none },
                    Config::Reader { "gadgetron_core_readers", "AcquisitionBucketReader", boost::none }
            };

            auto configs = default_readers;
            configs.insert(configs.end(), config.readers.begin(), config.readers.end());
            return load_readers(configs);
        }

        template<class CONFIG>
        std::vector<std::unique_ptr<Writer>> load_writers(CONFIG config) {
            return load_writers(config.writers);
        }

        template<class CONFIG>
        std::vector<std::unique_ptr<Writer>> load_default_and_additional_writers(CONFIG config) {

            static const std::vector<Config::Writer> default_writers{
                    Config::Writer { "gadgetron_mricore", "GadgetIsmrmrdAcquisitionMessageWriter" },
                    Config::Writer { "gadgetron_mricore", "GadgetIsmrmrdWaveformMessageWriter" },
                    Config::Writer { "gadgetron_core_writers", "ImageWriter" },
                    Config::Writer { "gadgetron_core_writers", "BufferWriter" },
                    Config::Writer { "gadgetron_core_writers", "IsmrmrdImageArrayWriter" },
                    Config::Writer { "gadgetron_core_writers", "AcquisitionBucketWriter" }
            };

            auto configs = default_writers;
            configs.insert(configs.begin(), config.writers.begin(), config.writers.end());
            return load_writers(configs);
        }

    private:
        boost::dll::shared_library load_library(const std::string &shared_library_name);

        const Core::StreamContext context;

        std::vector<boost::dll::shared_library> libraries = std::vector<boost::dll::shared_library>();
    };
}

