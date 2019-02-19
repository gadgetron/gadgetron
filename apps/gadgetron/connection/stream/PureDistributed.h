//
// Created by dchansen on 2/18/19.
//

#pragma once
#include "Processable.h"
#include "connection/Config.h"
#include "Reader.h"
#include "Writer.h"
#include <map>
#include "connection/Loader.h"

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
    };

}
