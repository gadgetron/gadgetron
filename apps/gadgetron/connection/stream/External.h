#pragma once

#include "connection/Config.h"
#include "connection/stream/Processable.h"
#include "connection/stream/Stream.h"

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

   };
}