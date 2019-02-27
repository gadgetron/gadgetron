#pragma once

#include <map>

#include "Processable.h"
#include "Reader.h"
#include "Writer.h"

#include "connection/Loader.h"
#include "connection/Config.h"

namespace Gadgetron::Server::Connection::Stream {

    class PureDistributed : public Processable {

    public:
        PureDistributed(const Config::PureDistributed& config, const Core::Context& context, Loader& loader);
        void process(Core::InputChannel input, Core::OutputChannel output, ErrorHandler& error_handler) override;
        const std::string& name() override;

    private:
        std::map<uint16_t, std::unique_ptr<Core::Reader>> readers;
        std::vector<std::unique_ptr<Core::Writer>> writers;
        Loader& loader;
        std::string remote_config;
        const Core::Context context;
        size_t nworkers=0;
    };

}
