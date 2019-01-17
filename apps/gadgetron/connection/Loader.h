#pragma once

#include "stream/Stream.h"
#include "Config.h"

#include "Context.h"
#include "Reader.h"
#include "Writer.h"

namespace Gadgetron::Server::Connection {


    class Loader {
        using Context = Core::Context;
        using GadgetProperties = Core::GadgetProperties;
        using Reader  = Gadgetron::Core::Reader;
        using Writer  = Gadgetron::Core::Writer;
        using Stream  = Stream::Stream;
    public:
        explicit Loader(const Context &);

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
            auto library = load_library(dll);
            return library.get_alias<FACTORY>(prefix + classname);
        }

    private:
        boost::dll::shared_library load_library(const std::string &shared_library_name);
        boost::filesystem::path make_library_path(const std::string &shared_library_name) const;

        const Context context;

        std::vector<boost::dll::shared_library> libraries = std::vector<boost::dll::shared_library>();
    };
}

