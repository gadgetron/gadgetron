#include "Parallel.h"

#include <map>
#include <memory>

#include "connection/Loader.h"

#include "Channel.h"
#include "Context.h"

namespace {
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Parallel;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Stream;

    using DecoratedBranch = Gadgetron::Server::Connection::Stream::Parallel::DecoratedBranch;
    using DecoratedMerge  = Gadgetron::Server::Connection::Stream::Parallel::DecoratedMerge;

    using branch_factory = Loader::generic_factory<Branch>;
    using merge_factory = Loader::generic_factory<Merge>;

    std::unique_ptr<DecoratedBranch> load_branch(const Config::Branch &conf, const Context &context, Loader &loader) {
        auto factory = loader.load_factory<branch_factory>("branch_factory_export_", conf.classname, conf.dll);
        return std::make_unique<DecoratedBranch>(factory(context, conf.properties), Config::name(conf));
    }

    std::unique_ptr<DecoratedMerge> load_merge(const Config::Merge &conf, const Context &context, Loader &loader) {
        auto factory = loader.load_factory<merge_factory>("merge_factory_export_", conf.classname, conf.dll);
        return std::make_unique<DecoratedMerge>(factory(context, conf.properties), Config::name(conf));
    }

    template<class KEY, class VALUE, class F>
    auto transform_map(std::map<KEY, VALUE>& input, F f) {
        using TRANSFORMED = std::remove_reference_t<decltype(f(input.begin()->second))>;
        auto result = std::map<KEY, TRANSFORMED>();

        for (auto &key_val : input) result.emplace(key_val.first, f(key_val.second));

        return result;
    }
}

namespace Gadgetron::Server::Connection::Stream {

    Parallel::Parallel(
            const Config::Parallel &config,
            const Core::Context &context,
            Loader &loader
    ) : branch(load_branch(config.branch, context, loader)), merge(load_merge(config.merge, context, loader)) {
        std::transform(
                config.streams.begin(), config.streams.end(),
                std::back_inserter(streams),
                [&](auto &stream_config) { return std::make_shared<Stream>(stream_config, context, loader); }
        );
    }

    void Parallel::process(
            InputChannel input,
            OutputChannel output,
            ErrorHandler &error_handler
    ) {

        std::vector<std::thread> threads;
        std::map<std::string, ChannelPair> input_channels;
        std::map<std::string, ChannelPair> output_channels;


        for (auto &stream : streams) {
            input_channels.emplace(stream->key, make_channel<MessageChannel>());
            output_channels.emplace(stream->key, make_channel<MessageChannel>());
        }

        threads.emplace_back(error_handler.run(
                [this](auto input, auto input_channels, auto output, auto &nested_handler) {
                    ErrorHandler handler{nested_handler, branch->key};
                    handler.handle([&](){branch->process(std::move(input), std::move(input_channels), std::move(output));});
                },
                std::move(input),
                transform_map(input_channels, [](auto& val) { return std::move(val.output); }),
                split(output),error_handler
        ));


        threads.emplace_back(error_handler.run(
                [this](auto output_channels, auto output, auto &nested_handler) {
                    ErrorHandler handler{nested_handler, branch->key};
                    handler.handle([&](){merge->process(std::move(output_channels), std::move(output));});
                },
                transform_map(output_channels, [](auto &val) { return std::move(val.input); }),
                std::move(output),error_handler
        ));

        for (auto &stream : streams) {
            threads.emplace_back(
                    Processable::process_async(stream, std::move(input_channels.at(stream->key).input),
                    std::move(output_channels.at(stream->key).output), error_handler)
            );
        }

        for (auto &thread : threads) { thread.join(); }
    }


    const std::string &Parallel::name() {
        static const std::string n = "Parallel";
        return n;
    }

    void Parallel::DecoratedBranch::process(
            InputChannel input,
            std::map<std::string, OutputChannel> output,
            OutputChannel bypass
    ) {
         branch->process(std::move(input), std::move(output), std::move(bypass));
    }

    Parallel::DecoratedBranch::DecoratedBranch(
            std::unique_ptr<Branch> branch,
            std::string key
    ) : branch(std::move(branch)), key(std::move(key)) {}

    void Parallel::DecoratedMerge::process(
            std::map<std::string, InputChannel> input,
            OutputChannel output
    ) {
         merge->process(std::move(input), std::move(output));
    }

    Parallel::DecoratedMerge::DecoratedMerge(
            std::unique_ptr<Merge> merge,
            std::string key
    ) : merge(std::move(merge)), key(std::move(key)) {}
}