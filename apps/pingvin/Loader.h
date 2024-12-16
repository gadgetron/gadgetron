#pragma once

#include <map>
#include <memory>

#include "Config.h"
#include "nodes/Stream.h"

#include "PropertyMixin.h"
#include "Context.h"

#include <boost/dll.hpp>

namespace Gadgetron::Main::Nodes {
    class Stream;
}

namespace Gadgetron::Main {

    class Loader {
        using Context = Core::Context;
        using Stream  = Gadgetron::Main::Nodes::Stream;

        using GadgetProperties = Core::GadgetProperties;
    public:
        explicit Loader(const Core::StreamContext &);

        std::unique_ptr<Stream> load(const Config::Stream &);

        template<class RESULT>
        using generic_factory = std::unique_ptr<RESULT>(
                const Context &,
                const std::string &,
                const GadgetProperties &
        );

        template<class FACTORY>
        FACTORY& load_factory(const std::string &prefix, const std::string &classname, const std::string &dll) {
            GINFO_STREAM("loading " << prefix << " - " << classname << " from the dll " << dll);
            auto library = load_library(dll);
            return library.get_alias<FACTORY>(prefix + classname);
        }

    private:
        boost::dll::shared_library load_library(const std::string &shared_library_name);

        const Core::StreamContext context;

        std::vector<boost::dll::shared_library> libraries = std::vector<boost::dll::shared_library>();
    };
}

