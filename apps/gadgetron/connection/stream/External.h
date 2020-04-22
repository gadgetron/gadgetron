#pragma once

#include <future>
#include <boost/asio.hpp>
#include <boost/process.hpp>

#include "connection/Config.h"
#include "connection/stream/Processable.h"
#include "connection/stream/Stream.h"

#include "connection/stream/common/ExternalChannel.h"
#include "connection/stream/common/Serialization.h"
#include "connection/stream/common/Configuration.h"

#include "parallel/Branch.h"
#include "parallel/Merge.h"

#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::Stream {

    class External : public Processable {
    public:
        using InputChannel = Core::GenericInputChannel;
        using OutputChannel = Core::OutputChannel;

        External(const Config::External &, const StreamContext &, Loader &);

        void process(
                InputChannel input,
                OutputChannel output,
                ErrorHandler &error_handler
        ) override;

        const std::string& name() override;

    private:
        std::shared_ptr<ExternalChannel> open_connection(Config::Execute, const Core::Context &);
        std::shared_ptr<ExternalChannel> open_connection(Config::Connect, const Core::Context &);
        std::shared_ptr<ExternalChannel> open_external_channel(const Config::External &, const Core::Context &);

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