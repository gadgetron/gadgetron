#include "Stream.h"

#include "Distributed.h"
#include "External.h"
#include "Parallel.h"
#include "ParallelProcess.h"
#include "PureDistributed.h"
#include "connection/core/Processable.h"

#include "connection/Loader.h"

#include "Node.h"

namespace {
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Nodes;
    using namespace std::string_literals;

    std::string print_action(const Config::Execute& execute){
        return "Execute block with name "s + execute.name + " of type " + execute.type;
    }

    std::string print_action(const Config::Connect& connect){
        return "Connect block on port "s + connect.port;
    }

    std::string print_action(const Config::Action& action){
        return visit([](auto& ac){return print_action(ac);}, action);
    }

    class NodeProcessable : public Processable {
    public:
        NodeProcessable(std::function<std::unique_ptr<Node>()> factory, std::string name) : factory(std::move(factory)), name_(std::move(name)) {}

        void process(GenericInputChannel input,
                OutputChannel output,
                ErrorHandler &
        ) override {
            auto node = factory();
            node->process(input, output);
        }

        const std::string& name() override {
            return name_;
        }

    private:
        std::function<std::unique_ptr<Node>()> factory;
        const std::string name_;
    };

    std::shared_ptr<Processable> load_node(const Config::Gadget &conf, const StreamContext &context, Loader &loader) {
        auto factory = loader.load_factory<Loader::generic_factory<Node>>("gadget_factory_export_", conf.classname,
                                                                          conf.dll);
        return std::make_shared<NodeProcessable>(
            [=]() {
                GDEBUG("Loading Gadget %s of class %s from %s\n", conf.name.c_str(), conf.classname.c_str(), conf.dll.c_str());
                return factory(context, conf.properties);
            },
            Config::name(conf)
        );
    }

    std::shared_ptr<Processable> load_node(const Config::Parallel &conf, const StreamContext &context, Loader &loader) {
        GDEBUG("Loading Parallel block\n");
        return std::make_shared<Nodes::Parallel>(conf, context, loader);
    }

    std::shared_ptr<Processable> load_node(const Config::External &conf, const StreamContext &context, Loader &loader) {
        GDEBUG("Loading External %s \n", print_action(conf.action).c_str());
        return std::make_shared<Nodes::External>(conf, context, loader);
    }

    std::shared_ptr<Processable> load_node(const Config::Distributed &conf, const StreamContext &context, Loader &loader) {
        GDEBUG("Loading Distributed block\n");
        return std::make_shared<Nodes::Distributed>(conf, context, loader);
    }

    std::shared_ptr<Processable> load_node(const Config::ParallelProcess& conf, const StreamContext& context, Loader& loader){
        GDEBUG("Loading ParalleProcess block\n");
        return std::make_shared<Nodes::ParallelProcess>(conf,context,loader);
    }

    std::shared_ptr<Processable> load_node(const Config::PureDistributed& conf, const StreamContext& context, Loader& loader){
        GDEBUG("Loading PureDistributed block\n");
        return std::make_shared<Nodes::PureDistributed>(conf,context,loader);
    }
}

namespace Gadgetron::Server::Connection::Nodes {

    Stream::Stream(const Config::Stream &config, const Core::StreamContext &context, Loader &loader) : key(config.key) {
        for (auto &node_config : config.nodes) {
            nodes.emplace_back(
                    Core::visit([&](auto n) { return load_node(n, context, loader); }, node_config)
            );
        }
    }

    void Stream::process(GenericInputChannel input,
            OutputChannel output,
            ErrorHandler &error_handler
    ) {
        if (empty()) return;

        std::vector<GenericInputChannel> input_channels{};
        input_channels.emplace_back(std::move(input));
        std::vector<OutputChannel> output_channels{};

        for (auto i = 0; i < nodes.size()-1; i++) {
            auto channel = make_channel<MessageChannel>();
            input_channels.emplace_back(std::move(channel.input));
            output_channels.emplace_back(std::move(channel.output));
        }

        output_channels.emplace_back(std::move(output));

        ErrorHandler nested_handler{error_handler, name()};

        std::vector<std::thread> threads(nodes.size());
        for (auto i = 0; i < nodes.size(); i++) {
            threads[i] = Processable::process_async(
                nodes[i],
                std::move(input_channels[i]),
                std::move(output_channels[i]),
                nested_handler
            );
        }

        for (auto &thread : threads) {
            thread.join();
        }
    }

    bool Stream::empty() const { return nodes.empty(); }
}

const std::string &Gadgetron::Server::Connection::Nodes::Stream::name() {
    static const std::string name = "Stream";
    return key.empty() ? name : key;
}
