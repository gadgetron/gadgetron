#ifndef GADGETRON_BUILDERS_H
#define GADGETRON_BUILDERS_H

#include <memory>

#include "Stream.h"
#include "Config.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Builders {

    class ReaderBuilder {
    public:
        ReaderBuilder(const Gadgetron::Server::Config &config);
        void process(std::function<void(std::unique_ptr<Gadgetron::Core::Reader>, boost::dll::shared_library)> on_reader);
    private:
        const Gadgetron::Server::Config &config;
    };

    class WriterBuilder {
    public:
        WriterBuilder(const Gadgetron::Server::Config &config);
        void process(std::function<void(std::unique_ptr<Gadgetron::Core::Writer>, boost::dll::shared_library)> on_writer);
    private:
        const Gadgetron::Server::Config &config;
    };

    class StreamBuilder {
    public:
        StreamBuilder(const Gadgetron::Server::Config &config, const Gadgetron::Core::Context &context);
        void process(std::function<void(std::unique_ptr<Gadgetron::Core::Node>, boost::dll::shared_library)> on_node);
    private:
        const Gadgetron::Server::Config &config;
        const Gadgetron::Core::Context &context;
    };
}

#endif //GADGETRON_BUILDERS_H
