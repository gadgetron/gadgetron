#ifndef GADGETRON_BUILDERS_H
#define GADGETRON_BUILDERS_H

#include <memory>

#include "Node.h"
#include "Reader.h"
#include "Writer.h"
#include "Config.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Builders {

    class Builder {
    public:
        explicit Builder(const Gadgetron::Server::Config &config, const Gadgetron::Core::Context::Paths &paths);
        boost::dll::shared_library load_library(const std::string &shared_library_name);
        boost::filesystem::path make_library_path(const std::string &shared_library_name);
    protected:
        const Gadgetron::Server::Config &config;
        const Gadgetron::Core::Context::Paths &paths;
    };

    class ReaderBuilder : public Builder {
        using reader_callback = std::function<void(uint16_t, std::unique_ptr<Gadgetron::Core::Reader>, boost::dll::shared_library)>;
    public:
        explicit ReaderBuilder(const Gadgetron::Server::Config &config, const Gadgetron::Core::Context::Paths &paths);
        void process(reader_callback on_reader);
    };

    class WriterBuilder : public Builder {
    public:
        explicit WriterBuilder(const Gadgetron::Server::Config &config, const Gadgetron::Core::Context::Paths &paths);
        void process(std::function<void(std::unique_ptr<Gadgetron::Core::Writer>, boost::dll::shared_library)> on_writer);
    };

    class StreamBuilder : public Builder {
        using node_callback = std::function<void(std::unique_ptr<Gadgetron::Core::Node>, boost::dll::shared_library)>;

    public:
        explicit StreamBuilder(
                const Gadgetron::Server::Config &config,
                const Gadgetron::Core::Context &context,
                std::shared_ptr<Gadgetron::Core::MessageChannel> input,
                std::shared_ptr<Gadgetron::Core::MessageChannel> output
        );
        void process(node_callback on_node);

        void load_node(const Gadgetron::Server::Config::Gadget &gadget_config, node_callback on_node);
        void load_node(const Gadgetron::Server::Config::Parallel &parallel_config, node_callback on_node);
        void load_node(const Gadgetron::Server::Config::Distributed &distributed_config, node_callback on_node);

    private:
        const Gadgetron::Core::Context &context;
        std::shared_ptr<Gadgetron::Core::MessageChannel> input, output;
    };
}

#endif //GADGETRON_BUILDERS_H
