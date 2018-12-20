#include "Loader.h"

#include <boost/range/combine.hpp>

#include "Connection.h"

#include "Node.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection;

namespace {

    class Stream : public Node {
    public:

        Stream(ErrorHandler &error_handler, std::vector<std::unique_ptr<Node>> nodes)
        : nodes(std::move(nodes)), error_handler(error_handler) {}

        void process(
                std::shared_ptr<InputChannel<Message>> in,
                std::shared_ptr<OutputChannel> out
        ) override {

            std::vector<std::shared_ptr<InputChannel<Message>>> input_channels(nodes.size());
            std::vector<std::shared_ptr<OutputChannel>> output_channels(nodes.size());

            input_channels.push_back(in);

            for (auto &node : nodes) {

                auto channel = std::make_shared<MessageChannel>();

                input_channels.push_back(channel);
                output_channels.push_back(channel);
            }

            output_channels.push_back(out);

            std::vector<std::thread> threads(nodes.size());
            for (auto i = 0; i < nodes.size(); i++) {
                threads[i] = start_node(*nodes[i], input_channels[i], output_channels[i]);
            }

            for (auto &thread : threads) {
                thread.join();
            }
        }

    private:

        std::thread start_node(
                Node &node,
                std::shared_ptr<InputChannel<Message>> input,
                std::shared_ptr<OutputChannel> output
        ) {
            std::function<void()> process = [&]() {
                node.process(input, output);
            };

            std::function<void()> handled_process = [&, process]() {
                // TODO: Node names.
                error_handler.handle("Unknown Node", process);
            };

            return std::thread(handled_process);
        }


        std::vector<std::unique_ptr<Node>> nodes;
        ErrorHandler &error_handler;
    };
}

namespace Gadgetron::Server::Connection {

    Loader::Loader(
            ErrorHandler &error_handler,
            Context context,
            Config config
    ) : error_handler(error_handler), context(std::move(context)), config(std::move(config)) {}

    boost::filesystem::path Loader::make_library_path(const std::string &shared_library_name) {
        return context.paths.gadgetron_home / "lib" / ("lib" + shared_library_name + ".so");
    }

    boost::dll::shared_library Loader::load_library(const std::string &shared_library_name) {

        std::lock_guard<std::mutex> guard(mutex);

        auto lib = boost::dll::shared_library(
                make_library_path(shared_library_name),
                boost::dll::load_mode::search_system_folders |
                boost::dll::load_mode::append_decorations
        );

        libraries.push_back(lib);
        return lib;
    }


    std::vector<std::pair<uint16_t, std::unique_ptr<Reader>>> Loader::readers() {

        std::vector<std::pair<uint16_t, std::unique_ptr<Reader>>> readers{};

        for (auto &reader_config : config.readers) {

            auto library = load_library(reader_config.dll);
            auto factory = library.get_alias<std::unique_ptr<Reader>(void)>(
                    "reader_factory_export_" + reader_config.classname);

            auto reader = factory();

            uint16_t slot = reader_config.slot.value_or(reader->slot());

            readers.emplace_back(std::make_pair(slot, std::move(reader)));
        }

        return std::move(readers);
    }

    std::vector<std::unique_ptr<Writer>> Loader::writers() {

        std::vector<std::unique_ptr<Writer>> writers{};

        for (auto &writer_config : config.writers) {

            auto library = load_library(writer_config.dll);
            auto factory = library.get_alias<std::unique_ptr<Writer>(void)>(
                    "writer_factory_export_" + writer_config.classname);

            auto writer= factory();

            writers.push_back(std::move(writer));
        }

        return std::move(writers);
    }

    std::unique_ptr<Node> Loader::stream() {

        std::vector<std::unique_ptr<Node>> nodes;
        for (auto &node_config : config.stream.nodes) {
            nodes.emplace_back(
                    boost::apply_visitor([&](auto n) { return load_node(n); }, node_config)
            );
        }

        return std::make_unique<Stream>(error_handler, std::move(nodes));
    }

    using gadget_factory = std::unique_ptr<Node>(
            const Context &,
            const std::unordered_map<std::string, std::string> &
    );

    std::unique_ptr<Node>
    Loader::load_node(const Config::Gadget &gadget_config) {
        auto library = load_library(gadget_config.dll);
        auto factory = library.get_alias<gadget_factory>("gadget_factory_export_" + gadget_config.classname);

        std::string name = gadget_config.name;
        if (name.empty()) name = gadget_config.dll;

        return factory(context, gadget_config.properties);
    }

    std::unique_ptr<Node>
    Loader::load_node(const Config::Parallel &parallel_config) {
        return std::unique_ptr<Node>();
    }

    std::unique_ptr<Node>
    Loader::load_node(const Config::Distributed &distributed_config) {
        return std::unique_ptr<Node>();
    }
}


