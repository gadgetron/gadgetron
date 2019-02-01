//
// Created by dchansen on 1/16/19.
//

#pragma once

#include "distributed/Distributor.h"
#include "Channel.h"
#include "connection/stream/Processable.h"
#include "Stream.h"
#include "connection/distributed/RemoteChannel.h"
#include <memory>
#include <map>
#include <vector>
#include <string>
#include "connection/distributed/CyclicIterator.h"

namespace Gadgetron::Server::Connection::Stream {
    class Distributed : public Processable {

    public:
        using RemoteChannel = Gadgetron::Server::Distributed::RemoteChannel;

        using Worker = Core::variant<Server::Distributed::Address, Server::Distributed::Local>;
        using Address = Gadgetron::Server::Distributed::Address;
        using Local = Gadgetron::Server::Distributed::Local;

        Distributed(const Config::Distributed &, const Core::Context &, Loader &);

        const std::string &name() override;


        void process(Core::InputChannel input, Core::OutputChannel output,
                     ErrorHandler &error_handler) override;

    private:

        std::unique_ptr<Core::Distributed::Distributor> load_distributor(const Config::Distributor &);


        std::unique_ptr<Core::Distributed::Distributor> distributor;
        Gadgetron::Server::Connection::Loader &loader;
        const Core::Context context;
        const Config::Distributed config;
    };
}



