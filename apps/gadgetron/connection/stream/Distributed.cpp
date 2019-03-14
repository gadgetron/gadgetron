#include <algorithm>

#include "Distributed.h"

#include "connection/stream/distributed/Discovery.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Server::Connection;

}

namespace Gadgetron::Server::Connection::Stream {

    void Distributed::process(
            Core::InputChannel input,
            Core::OutputChannel output,
            ErrorHandler& error_handler
    ) {
//        auto channel_creator = WorkerChannelCreator{
//            config,
//            context,
//            loader,
//            Core::split(output),
//            error_handler
//        };
//
//
//        distributor->process(std::move(input), channel_creator, std::move(output));
    }

    Distributed::Distributed(
            const Config::Distributed& distributed_config,
            const Core::Context& context,
            Loader& loader
    ) : context{ context },
        loader{ loader },
        config(distributed_config) {
        distributor = load_distributor(distributed_config.distributor);
    }

    std::unique_ptr<Core::Distributed::Distributor> Distributed::load_distributor(
            const Config::Distributor& conf
    ) {
        auto factory = loader.load_factory<Loader::generic_factory<Core::Distributed::Distributor>>(
                "distributor_factory_export_", conf.classname, conf.dll);
        return factory(context, conf.properties);
    }

    const std::string& Distributed::name() {
        static const std::string n = "Distributed";
        return n;
    }
}

