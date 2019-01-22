//
// Created by dchansen on 1/16/19.
//

#pragma once

#include "distributed/Distributor.h"
#include "Channel.h"
#include "connection/stream/Processable.h"
#include "Stream.h"
#include "connection/distributed/RemoteChannel.h"

namespace Gadgetron::Server::Connection::Stream {
    class Distributed : public Processable {

        using RemoteChannel = Gadgetron::Server::Distributed::RemoteChannel;
        Distributed(const Config::Distributed &, const Core::Context &, Loader &);


    public:
        void process(std::shared_ptr<Core::Channel> input, std::shared_ptr<Core::Channel> output,
                     ErrorHandler &error_handler) override;

    private:

        std::shared_ptr<RemoteChannel> create_remote_channel();

        std::unique_ptr<Core::Distributed::Distributor> load_distributor(const Config::Distributor&);

        std::unique_ptr<Core::Distributed::Distributor> distributor;
        const Config config;
        const Core::Context context;
        std::map<uint16_t,std::unique_ptr<Core::Reader>> readers;
        std::vector<std::unique_ptr<Core::Writer>> writers;
        Gadgetron::Server::Connection::Loader& loader;
        const std::vector<Gadgetron::Server::Distributed::Address> workers;

    };
}



