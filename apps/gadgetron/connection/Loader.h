#pragma once

#include "Core.h"
#include "Config.h"

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
        using Node    = Gadgetron::Core::Node;

    public:
        Loader(ErrorHandler &error_handler, Context context, Config config);

        std::vector<std::pair<std::uint16_t, std::unique_ptr<Reader>>> readers();
        std::vector<std::unique_ptr<Writer>> writers();
        std::unique_ptr<NodeHandler> stream();

    private:
        boost::filesystem::path make_library_path(const std::string &shared_library_name);
        boost::dll::shared_library load_library(const std::string &shared_library_name);

        std::unique_ptr<NodeHandler> load_node(const Config::Gadget& gadget_config);
        std::unique_ptr<NodeHandler> load_node(const Config::Parallel& parallel_config);
        std::unique_ptr<NodeHandler> load_node(const Config::Distributed& distributed_config);

        ErrorHandler &error_handler;

        const Context context;
        const Config config;

        std::mutex mutex;
        std::vector<boost::dll::shared_library> libraries;
    };

}

