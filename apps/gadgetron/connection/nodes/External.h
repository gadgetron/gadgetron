#pragma once

#include <future>
#include <boost/asio.hpp>
#include <boost/process.hpp>

#include "connection/config/Config.h"

#include "Stream.h"
#include "connection/core/Processable.h"

#include "common/ExternalChannel.h"
#include "common/Serialization.h"
#include "common/Configuration.h"

#include "parallel/Branch.h"
#include "parallel/Merge.h"

#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::Nodes {

    class External : public Processable {
    public:
        using InputChannel = Core::GenericInputChannel;
        using OutputChannel = Core::OutputChannel;

        External(const Config::External &, const Core::StreamContext &, Loader &);

        void process(
                InputChannel input,
                OutputChannel output,
                ErrorHandler &error_handler
        ) override;

        const std::string& name() override;

    private:
        std::shared_ptr<ExternalChannel> open_connection(Config::Execute, const Core::StreamContext &);
        std::shared_ptr<ExternalChannel> open_connection(Config::Connect, const Core::StreamContext &);
        std::shared_ptr<ExternalChannel> open_external_channel(const Config::External &, const Core::StreamContext &);

        void monitor_child(std::shared_ptr<boost::process::child>, std::shared_ptr<boost::asio::ip::tcp::acceptor>);

        std::future<std::shared_ptr<ExternalChannel>> channel;
        std::shared_ptr<Serialization> serialization;
        std::shared_ptr<Configuration> configuration;

        boost::asio::io_service io_service;

        struct {
            std::future<void> child;
        } monitors;

   };
}