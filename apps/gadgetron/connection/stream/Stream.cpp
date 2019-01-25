#include "Stream.h"

#include "connection/stream/Processable.h"
#include "connection/stream/Parallel.h"
#include "connection/stream/External.h"
#include "connection/stream/Distributed.h"
#include "connection/CloseGuard.h"
#include "connection/Loader.h"

#include "Node.h"

namespace {
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Stream;

    class NodeProcessable : public Processable {
    public:
        NodeProcessable(std::unique_ptr<Node> node, std::string name) : node(std::move(node)), name(std::move(name)) {}

        void process(
                std::shared_ptr<Channel> input,
                std::shared_ptr<Channel> output,
                ErrorHandler &error_handler
        ) override {
            CloseGuard closer{input, output};
            error_handler.handle(name, [&]() { node->process(*input, *output); });
        }

    private:
        std::unique_ptr<Node> node;
        const std::string name;
    };

    std::unique_ptr<Processable> load_node(const Config::Gadget &conf, const Context &context, Loader &loader) {
        auto factory = loader.load_factory<Loader::generic_factory<Node>>("gadget_factory_export_", conf.classname, conf.dll);
        return std::make_unique<NodeProcessable>(factory(context, conf.properties), Config::name(conf));
    }

    std::unique_ptr<Processable> load_node(const Config::Parallel &conf, const Context &context, Loader &loader) {
        return std::make_unique<Gadgetron::Server::Connection::Stream::Parallel>(conf, context, loader);
    }

    std::unique_ptr<Processable> load_node(const Config::Distributed &conf, const Context &context, Loader &loader) {
        return std::make_unique<Gadgetron::Server::Connection::Stream::Distributed>(conf,context,loader);
    }
}

namespace Gadgetron::Server::Connection::Stream {

    Stream::Stream(const Config::Stream &config, const Core::Context &context, Loader &loader) : key(config.key) {
        for (auto &node_config : config.nodes) {
            nodes.emplace_back(
                    boost::apply_visitor([&](auto n) { return load_node(n, context, loader); }, node_config)
            );
        }
    }

    void Stream::process(
            std::shared_ptr<Channel> input,
            std::shared_ptr<Channel> output,
            ErrorHandler &error_handler
    ) {
        std::vector<std::shared_ptr<Channel>> input_channels{};
        std::vector<std::shared_ptr<Channel>> output_channels{};

        input_channels.push_back(input);

        for (auto i = 0; i < (nodes.size() - 1); i++) {

            auto channel = std::make_shared<MessageChannel>();

            input_channels.push_back(channel);
            output_channels.push_back(channel);
        }

        output_channels.push_back(output);

        DecoratedErrorHandler nested_handler{error_handler, key};

        std::vector<std::thread> threads(nodes.size());
        for (auto i = 0; i < nodes.size(); i++) {
            threads[i] = error_handler.run(
                    key,
                    [&, i](auto in, auto out) {
                        nodes[i]->process(in, out, nested_handler);
                    },
                    input_channels[i],
                    output_channels[i]
            );
        }

        for (auto &thread : threads) {
            thread.join();
        }
    }
}
