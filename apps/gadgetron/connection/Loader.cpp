#include "Loader.h"

#include <map>
#include <memory>
#include <boost/range/combine.hpp>

#include "Connection.h"
#include "NodeHandler.h"

#include "parallel/Branch.h"
#include "parallel/Merge.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Parallel;
using namespace Gadgetron::Server::Connection;

namespace {

    class DecoratedErrorHandler : public ErrorHandler {
    public:
        ErrorHandler &handler;
        std::string decorator;

        DecoratedErrorHandler(ErrorHandler &handler, std::string decorator)
        : handler(handler), decorator(std::move(decorator)) {}

        void handle(const std::string &location, std::function<void()> function) override {
            handler.handle(decorator + "/" + location, function);
        }
    };

    class NodeHousekeeping : public NodeHandler {
    public:
        NodeHousekeeping(
                std::unique_ptr<Node> node,
                ErrorHandler &error_handler,
                std::string location
        ) : node(std::move(node)), location(std::move(location)), error_handler(error_handler) {};

        void process(
                std::shared_ptr<Channel> in,
                std::shared_ptr<Channel> out
        ) override {
            error_handler.handle(location, [&]() {
                node->process(*in, *out);
            });
            in->close();
            out->close();
        }

    private:
        std::unique_ptr<Node> node;
        const std::string location;
        ErrorHandler &error_handler;
    };

    class BranchHousekeeping : public BranchHandler {
    public:
        BranchHousekeeping(
                std::unique_ptr<Branch> branch,
                ErrorHandler &error_handler
        ) {};
    };

    class MergeHousekeeping : public MergeHandler {
    public:
        MergeHousekeeping(
                std::unique_ptr<Merge> merge,
                ErrorHandler &error_handler
        ) {};


    };

    class Stream : public NodeHandler {
    public:
        Stream(std::vector<std::unique_ptr<NodeHandler>> nodes) : nodes(std::move(nodes)) {}

        void process(
                std::shared_ptr<Channel> in,
                std::shared_ptr<Channel> out
        ) override {

            std::vector<std::shared_ptr<Channel>> input_channels{};
            std::vector<std::shared_ptr<Channel>> output_channels{};

            input_channels.push_back(in);

            for (auto i = 0; i < (nodes.size() - 1); i++) {

                auto channel = std::make_shared<MessageChannel>();

                input_channels.push_back(channel);
                output_channels.push_back(channel);
            }

            output_channels.push_back(out);

            std::vector<std::thread> threads(nodes.size());
            for (auto i = 0; i < nodes.size(); i++) {
                threads[i] = std::thread(
                        [&, i](auto in, auto out) {
                            nodes[i]->process(in, out);
                        },
                        input_channels[i],
                        output_channels[i]
                );
            }

            for (auto &thread : threads) {
                thread.join();
            }
        }

    private:
        std::vector<std::unique_ptr<NodeHandler>> nodes;
    };

    class ParallelNode : public NodeHandler {
    public:

        class ParallelStream {
        public:
            struct {
                std::shared_ptr<MessageChannel> input, output;
            } channels;

            std::string key;
            std::unique_ptr<NodeHandler> stream;
            std::unique_ptr<ErrorHandler> error_handler;

            ParallelStream(
                    std::string key,
                    std::unique_ptr<NodeHandler> stream,
                    ErrorHandler &error_handler
            ) : key(std::move(key)), stream(std::move(stream)),
                error_handler(std::make_unique<DecoratedErrorHandler>(error_handler, key)),
                channels{ std::make_shared<MessageChannel>(), std::make_shared<MessageChannel>() } {}

            std::thread start_process_thread() {
                return error_handler->run(
                    key,
                    [=]() { stream->process(channels.input, channels.output); }
                );
            }
        };


        ParallelNode(
                std::unique_ptr<Branch> branch,
                std::unique_ptr<Merge>  merge,
                std::map<std::string, std::unique_ptr<NodeHandler>> streams,
                ErrorHandler &error_handler
        ) : branch(std::move(branch)), merge(std::move(merge)), error_handler(error_handler) {

            for (auto &pair : streams) {
                this->streams.emplace_back(std::make_unique<ParallelStream>(
                        pair.first,
                        std::move(pair.second),
                        error_handler
                ));
            }
        }

        void process(
                std::shared_ptr<Channel> in,
                std::shared_ptr<Channel> out
        ) override {
            std::vector<std::thread> threads;
            std::map<std::string, std::shared_ptr<Channel>> input_channels, output_channels;

            std::transform(streams.begin(), streams.end(), std::inserter(input_channels, input_channels.begin()),
                    [](auto &stream) { return std::make_pair(stream->key, stream->channels.input); }
            );

            std::transform(streams.begin(), streams.end(), std::inserter(output_channels, output_channels.begin()),
                    [](auto &stream) { return std::make_pair(stream->key, stream->channels.output); }
            );

            threads.emplace_back(error_handler.run(
                    "Branch",
                    [&]() { branch->process(in, input_channels, out); }
            ));

            threads.emplace_back(error_handler.run(
                    "Merge",
                    [&]() { merge->process(output_channels, out); }
            ));

            for (auto &stream : streams) {
                threads.emplace_back(stream->start_process_thread());
            }

            for(auto &thread : threads) { thread.join(); }
        }

    private:
        std::unique_ptr<Branch> branch;
        std::unique_ptr<Merge>  merge;
        std::vector<std::unique_ptr<ParallelStream>> streams;

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
        return context.paths.gadgetron_home / "lib" / shared_library_name;
    }

    boost::dll::shared_library Loader::load_library(const std::string &shared_library_name) {

        auto lib = boost::dll::shared_library(
                make_library_path(shared_library_name),
                boost::dll::load_mode::append_decorations
        );

        libraries.push_back(lib);
        return lib;
    }

    std::vector<std::pair<uint16_t, std::unique_ptr<Reader>>> Loader::readers() {

        std::vector<std::pair<uint16_t, std::unique_ptr<Reader>>> readers{};

        for (auto &reader_config : config.readers) {
            auto reader = load_reader(reader_config);
            uint16_t slot = reader_config.slot.value_or(reader->slot());
            readers.emplace_back(std::make_pair(slot, std::move(reader)));
        }

        return std::move(readers);
    }

    std::vector<std::unique_ptr<Writer>> Loader::writers() {

        std::vector<std::unique_ptr<Writer>> writers{};

        for (auto &writer_config : config.writers) {
            writers.emplace_back(load_writer(writer_config));
        }

        return std::move(writers);
    }

    std::unique_ptr<NodeHandler> Loader::stream() {
        return load_stream(config.stream);
    }

    std::unique_ptr<Reader> Loader::load_reader(const Config::Reader &reader_config) {
        auto library = load_library(reader_config.dll);
        auto factory = library.get_alias<std::unique_ptr<Reader>(void)>(
                "reader_factory_export_" + reader_config.classname);
        return factory();
    }

    std::unique_ptr<Writer> Loader::load_writer(const Config::Writer &writer_config) {
        auto library = load_library(writer_config.dll);
        auto factory = library.get_alias<std::unique_ptr<Writer>(void)>(
                "writer_factory_export_" + writer_config.classname);
        return factory();
    }

    std::unique_ptr<BranchHandler> Loader::load_branch(const Config::Branch &branch_config) {

        auto library = load_library(branch_config.dll);
        auto factory = library.get_alias<branch_factory>(
                "branch_factory_export_" + branch_config.classname);

        return std::make_unique<BranchHousekeeping>(
                factory(context, branch_config.properties)
        );
    }

    std::unique_ptr<MergeHandler> Loader::load_merge(const Config::Merge &merge_config) {

        auto library = load_library(merge_config.dll);
        auto factory = library.get_alias<merge_factory>(
                "merge_factory_export_" + merge_config.classname);

        return std::make_unique<MergeHousekeeping>(
                factory(context, merge_config.properties)
        );
    }

    std::unique_ptr<NodeHandler> Loader::load_stream(const Config::Stream &stream_config) {
        std::vector<std::unique_ptr<NodeHandler>> nodes;

        for (auto &node_config : stream_config.nodes) {
            nodes.emplace_back(
                    boost::apply_visitor([&](auto n) { return load_node(n); }, node_config)
            );
        }

        return std::make_unique<Stream>(std::move(nodes));
    }

    std::unique_ptr<NodeHandler> Loader::load_node(const Config::Gadget &gadget_config) {
        auto library = load_library(gadget_config.dll);
        auto factory = library.get_alias<node_factory>("gadget_factory_export_" + gadget_config.classname);

        std::string name = gadget_config.name;
        if (name.empty()) name = gadget_config.classname;

        return std::make_unique<NodeHousekeeping>(
                factory(context, gadget_config.properties),
                error_handler,
                name
        );
    }

    std::unique_ptr<NodeHandler> Loader::load_node(const Config::Parallel &parallel_config) {

        unsigned int stream_index = 0;

        std::map<std::string, std::unique_ptr<NodeHandler>> streams;
        for (auto &stream_config : parallel_config.streams) {

            std::string name = stream_config.name;
            if (name.empty()) name = "unnamed-stream-" + std::to_string(stream_index++);

            streams[name] = load_stream(stream_config);
        }

        return std::make_unique<ParallelNode>(
                load_branch(parallel_config.branch),
                load_merge(parallel_config.merge),
                std::move(streams),
                error_handler
        );
    }

    std::unique_ptr<NodeHandler> Loader::load_node(const Config::Distributed &distributed_config) {
        return std::unique_ptr<NodeHandler>();
    }
}
