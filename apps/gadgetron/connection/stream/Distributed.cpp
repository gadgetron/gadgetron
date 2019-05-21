
#include <list>
#include <algorithm>

#include "Distributed.h"

#include "connection/stream/distributed/Discovery.h"
#include "connection/stream/distributed/Worker.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Distributed;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Stream;

    class ChannelWrapper {
    public:
        ChannelWrapper(
                Address peer,
                std::shared_ptr<Serialization> serialization,
                std::shared_ptr<Configuration> configuration,
                InputChannel input,
                OutputChannel output
        );

    private:
        Address peer;

        std::shared_ptr<Serialization> serialization;
        std::shared_ptr<Configuration> configuration;

        std::unique_ptr<std::iostream> stream;

        InputChannel input;
        OutputChannel output;
    };

    ChannelWrapper::ChannelWrapper(
            Address peer,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration,
            InputChannel input,
            OutputChannel output
    ) : peer(std::move(peer)),
        serialization(std::move(serialization)),
        configuration(std::move(configuration)),
        input(std::move(input)),
        output(std::move(output)) {

        stream = connect(peer, configuration);
    }


    class ChannelCreatorImpl : public ChannelCreator {
    public:
        OutputChannel create() override;

        ChannelCreatorImpl(
                std::shared_ptr<Serialization> serialization,
                std::shared_ptr<Configuration> configuration,
                OutputChannel output_channel,
                ErrorHandler& error_handler
        );

    private:
        Address next_peer();

        OutputChannel output;

        std::shared_ptr<Serialization> serialization;
        std::shared_ptr<Configuration> configuration;

        std::list<Address> peers;

        ErrorHandler error_handler;
    };

    ChannelCreatorImpl::ChannelCreatorImpl(
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration,
            OutputChannel output_channel,
            ErrorHandler &error_handler
    ) : serialization(std::move(serialization)),
        configuration(std::move(configuration)),
        output(std::move(output_channel)),
        error_handler(error_handler, "Distributed") {

        auto ps = discover_peers();
        std::copy(ps.begin(), ps.end(), std::back_inserter(peers));
    }

    OutputChannel ChannelCreatorImpl::create() {

        auto pair = Core::make_channel<MessageChannel>();

        auto channel = ChannelWrapper(
                next_peer(),
                serialization,
                configuration,
                std::move(pair.input),
                Core::split(output)
        );

        return std::move(pair.output);
    }

    Address ChannelCreatorImpl::next_peer() {
        auto peer = peers.front(); peers.pop_front();
        peers.push_back(peer);
        return peer;
    }
}

namespace {
    std::unique_ptr<Core::Distributed::Distributor> load_distributor(
            Loader &loader,
            const Core::Context& context,
            const Config::Distributor& conf
    ) {
        auto factory = loader.load_factory<Loader::generic_factory<Core::Distributed::Distributor>>(
                "distributor_factory_export_", conf.classname, conf.dll);
        return factory(context, conf.properties);
    }
}

namespace Gadgetron::Server::Connection::Stream {

    void Distributed::process(
            Core::InputChannel input,
            Core::OutputChannel output,
            ErrorHandler& error_handler
    ) {
        auto channel_creator = ChannelCreatorImpl {
            serialization,
            configuration,
            Core::split(output),
            error_handler
        };

        distributor->process(std::move(input), channel_creator, std::move(output));
    }

    Distributed::Distributed(
            const Config::Distributed& config,
            const Core::Context& context,
            Loader& loader
    ) : serialization(std::make_shared<Serialization>(
                loader.load_readers(config),
                loader.load_writers(config)
        )),
        configuration(std::make_shared<Configuration>(
                context,
                config
        )),
        distributor(load_distributor(loader, context, config.distributor)) {}

    const std::string& Distributed::name() {
        static const std::string n = "Distributed";
        return n;
    }
}

