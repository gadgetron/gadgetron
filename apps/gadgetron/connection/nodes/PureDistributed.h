#pragma once

#include <map>
#include <list>
#include <thread>

#include "Reader.h"
#include "Writer.h"
#include "connection/core/Processable.h"

#include "connection/Loader.h"
#include "connection/config/Config.h"

#include "common/Serialization.h"
#include "common/Configuration.h"

#include "distributed/Worker.h"
#include "distributed/Pool.h"

namespace Gadgetron::Server::Connection::Nodes {

    class PureDistributed : public Processable {

    public:
        PureDistributed(
                const Config::PureDistributed &,
                const Core::StreamContext &,
                Loader &
        );

        void process(
                Core::GenericInputChannel,
                Core::OutputChannel,
                ErrorHandler &
        ) override;

        const std::string& name() override;

    private:
        using Job = std::future<Core::Message>;
        using Queue = Core::MPMCChannel<Job>;

        void process_outbound(Core::GenericInputChannel, std::shared_ptr<Queue>);
        void process_inbound(Core::OutputChannel, std::shared_ptr<Queue>);

        std::shared_ptr<Serialization> serialization;
        std::shared_ptr<Configuration> configuration;

        std::list<std::future<std::unique_ptr<Worker>>> pending_workers;
    };
}
