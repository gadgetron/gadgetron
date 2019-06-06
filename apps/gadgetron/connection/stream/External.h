#pragma once

#include <future>
#include <boost/asio.hpp>

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
        using InputChannel = Core::InputChannel;
        using OutputChannel = Core::OutputChannel;

        External(const Config::External &, const Core::Context &, Loader &);

        void process(
                InputChannel input,
                OutputChannel output,
                ErrorHandler &error_handler
        ) override;

        const std::string& name() override;

    private:
        std::unique_ptr<std::iostream> open_connection(Config::Execute, const Core::Context &);
        std::unique_ptr<std::iostream> open_connection(Config::Connect, const Core::Context &);
        std::shared_ptr<ExternalChannel> open_external_channel(const Config::External &, const Core::Context &);

        std::future<std::shared_ptr<ExternalChannel>> channel;
        std::shared_ptr<Serialization> serialization;
        std::shared_ptr<Configuration> configuration;

        boost::asio::io_service io_service;
   };
}