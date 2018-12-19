#ifndef GADGETRON_BUILDERS_H
#define GADGETRON_BUILDERS_H

#include <memory>

#include "Node.h"
#include "Reader.h"
#include "Writer.h"
#include "Config.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {

    class Builder {
    public:
        Builder(const Gadgetron::Core::Context::Paths& context_in);

        std::vector<std::unique_ptr<Gadgetron::Core::Writer>> build_writers(const std::vector<Config::Writer>& writers);
        std::vector<std::pair<uint16_t,std::unique_ptr<Gadgetron::Core::Reader>>> build_readers(const std::vector<Config::Reader>& readers);
        std::unique_ptr<Core::Node> build_stream(const Config::Stream& stream, const Core::Context& context);

    private:

        boost::dll::shared_library load_library(const std::string &shared_library_name);
        std::unique_ptr<Gadgetron::Core::Node> load_node(const Config::Gadget& gadget_config, const Core::Context& context);
        std::unique_ptr<Gadgetron::Core::Node> load_node(const Config::Parallel& parallel_config, const Core::Context& );
        std::unique_ptr<Gadgetron::Core::Node> load_node(const Config::Distributed& distributed_config, const Core::Context& );

        std::vector<boost::dll::shared_library> libraries;
        const Gadgetron::Core::Context::Paths paths;
    };




}

#endif //GADGETRON_BUILDERS_H
