#pragma once

#include <map>
#include <vector>
#include <string>
#include <memory>

#include "distributed/Distributor.h"

#include "Processable.h"

#include "Channel.h"
#include "Stream.h"

namespace Gadgetron::Server::Connection::Stream {

    class Distributed : public Processable {
    public:
        Distributed(
                const Config::Distributed &,
                const Core::Context &,
                Loader &
        );

        void process(
                Core::InputChannel input,
                Core::OutputChannel output,
                ErrorHandler &error_handler
        ) override;

        const std::string &name() override;

    private:
        std::unique_ptr<Core::Distributed::Distributor> load_distributor(const Config::Distributor &);
        std::unique_ptr<Core::Distributed::Distributor> distributor;

        const Core::Context context;
        const Config::Distributed config;

        Gadgetron::Server::Connection::Loader &loader;
    };
}



