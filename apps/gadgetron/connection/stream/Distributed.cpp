
#include <list>
#include <algorithm>

#include "Distributed.h"

#include "connection/stream/common/Closer.h"
#include "connection/stream/common/Discovery.h"
#include "connection/stream/common/ExternalChannel.h"
#include "io/iostream_operators.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Distributed;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Stream;

    class ChannelWrapper {
    public:
        ChannelWrapper(
                const Address &peer,
                std::shared_ptr<Serialization> serialization,
                std::shared_ptr<Configuration> configuration
        );

        void process_input(GenericInputChannel input);
        void process_output(OutputChannel output);

    private:
        std::shared_ptr<ExternalChannel> external;
    };

    ChannelWrapper::ChannelWrapper(
            const Address &peer,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration
    ) {
        GINFO_STREAM("Connecting to peer: " << peer);
        external = std::make_shared<ExternalChannel>(
                connect(peer, configuration),
                std::move(serialization),
                std::move(configuration)
        );
    }

    void ChannelWrapper::process_input(GenericInputChannel input) {
        auto closer = make_closer(external);
        for (auto message : input) {
            external->push_message(std::move(message));
        }
    }

    void ChannelWrapper::process_output(OutputChannel output) {
        while(true) {
            output.push_message(external->pop());
            GDEBUG_STREAM("Pushed message to distributed output.");
        }
    }

    class ChannelCreatorImpl : public ChannelCreator {
    public:
        OutputChannel create() override;
        void join();

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
        std::list<std::thread> threads;

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

        auto channel = std::make_shared<ChannelWrapper>(
                next_peer(),
                serialization,
                configuration
        );

        threads.push_back(error_handler.run(
                [=](auto in) { channel->process_input(std::move(in)); },
                std::move(pair.input)
        ));
        threads.push_back(error_handler.run(
                [=](auto out) { channel->process_output(std::move(out)); },
                Core::split(output)
        ));

        return std::move(pair.output);
    }

    void ChannelCreatorImpl::join() {
        for (auto &thread : threads) thread.join();
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
            Core::GenericInputChannel input,
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
        channel_creator.join();
    }

    Distributed::Distributed(
            const Config::Distributed& config,
            const Core::StreamContext& context,
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

