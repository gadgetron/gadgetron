#pragma once

#include "Core.h"
#include "Config.h"

#include "parallel/Merge.h"
#include "parallel/Branch.h"
#include "Context.h"
#include "Reader.h"
#include "Writer.h"
#include "Node.h"
#include "NodeHandler.h"


namespace Gadgetron::Server::Connection {

    class Loader {
        using Context = Gadgetron::Core::Context;
        using Reader  = Gadgetron::Core::Reader;
        using Writer  = Gadgetron::Core::Writer;
        using GadgetProperties = Gadgetron::Core::GadgetProperties;

        template<class RESULT>
        using generic_factory = std::unique_ptr<RESULT>(
                const Context &,
                const GadgetProperties &
        );
        using node_factory = generic_factory<Core::Node>;
        using branch_factory = generic_factory<Core::Parallel::Branch>;
        using merge_factory = generic_factory<Core::Parallel::Merge>;

    public:
        Loader(ErrorHandler &error_handler, Context context, Config config);

        std::vector<std::pair<std::uint16_t, std::unique_ptr<Reader>>> readers();
        std::vector<std::unique_ptr<Writer>> writers();
        std::unique_ptr<NodeHandler> stream();

    private:
        boost::filesystem::path make_library_path(const std::string &shared_library_name);
        boost::dll::shared_library load_library(const std::string &shared_library_name);

        std::unique_ptr<Reader> load_reader(const Config::Reader &);
        std::unique_ptr<Writer> load_writer(const Config::Writer &);
        std::unique_ptr<BranchHandler> load_branch(const Config::Branch &);
        std::unique_ptr<MergeHandler>  load_merge(const Config::Merge &);

        std::unique_ptr<NodeHandler> load_stream(const Config::Stream &);
        std::unique_ptr<NodeHandler> load_node(const Config::Gadget &);
        std::unique_ptr<NodeHandler> load_node(const Config::Parallel &);
        std::unique_ptr<NodeHandler> load_node(const Config::Distributed &);

        const Context context;
        const Config config;

        ErrorHandler &error_handler;

        std::vector<boost::dll::shared_library> libraries = std::vector<boost::dll::shared_library>();
    };
}

