#pragma once

#include <map>
#include <vector>
#include <string>
#include <memory>

#include "common/Configuration.h"
#include "common/Serialization.h"
#include "distributed/Distributor.h"

#include "Processable.h"

#include "Channel.h"
#include "Stream.h"

namespace Gadgetron::Server::Connection::Stream {

    class Distributed : public Processable {
    public:
        Distributed(
                const Config::Distributed &,
                const Core::StreamContext &,
                Loader &
        );

        void process(
                Core::GenericInputChannel input,
                Core::OutputChannel output,
                ErrorHandler &error_handler
        ) override;

        const std::string &name() override;

    private:
        std::unique_ptr<Core::Distributed::Distributor> distributor;

        const std::shared_ptr<Serialization> serialization;
        const std::shared_ptr<Configuration> configuration;
    };
}



