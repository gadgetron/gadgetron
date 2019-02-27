#pragma once

#include "PureStream.h"
#include "connection/stream/Processable.h"

namespace Gadgetron::Server::Connection::Stream {
    class ParallelProcess : public Processable {

    public:
        ParallelProcess(const Config::ParallelProcess& conf, const Core::Context& context, Loader& loader);
        void process(Core::InputChannel input, Core::OutputChannel output, ErrorHandler& error_handler) override;
        const std::string& name() override;
    private:
        const size_t workers;
        const PureStream pureStream;
    };
}
