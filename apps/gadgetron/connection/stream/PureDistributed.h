#pragma once

#include <map>
#include <list>

#include "Processable.h"
#include "Reader.h"
#include "Writer.h"

#include "connection/Loader.h"
#include "connection/Config.h"

#include "common/Serialization.h"
#include "common/Configuration.h"

#include "distributed/Worker.h"
#include "distributed/Pool.h"

namespace Gadgetron::Server::Connection::Stream {

    class PureDistributed : public Processable {

    public:
        PureDistributed(
                const Config::PureDistributed &,
                const Core::Context &,
                Loader &
        );

        void process(
                Core::InputChannel,
                Core::OutputChannel,
                ErrorHandler &
        ) override;

        const std::string& name() override;

    private:
        struct Job;
        using Queue = Core::MPMCChannel<Job>;

        Job send_message_to_worker(Core::Message message, std::shared_ptr<Worker> worker);
        Core::Message get_message_from_worker(Job job, size_t retries = 3);

        void process_outbound(Core::InputChannel, Queue &);
        void process_inbound(Core::OutputChannel, Queue &);

        void trigger_node_discovery();

        const std::shared_ptr<Serialization> serialization;
        const std::shared_ptr<Configuration> configuration;

        std::shared_ptr<Pool<Worker>> workers;
    };
}
