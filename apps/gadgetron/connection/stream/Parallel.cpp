#include "Parallel.h"

#include <map>
#include <memory>

#include "connection/Loader.h"
#include "connection/CloseGuard.h"

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
                [&](auto &stream_config) { return std::make_unique<Stream>(stream_config, context, loader) ; }
        );
    }

    void Parallel::process(
            std::shared_ptr<Channel> input,
            std::shared_ptr<Channel> output,
            ErrorHandler &error_handler
    ) {

        DecoratedErrorHandler nested_handler{error_handler, "parallel"};

        std::vector<std::thread> threads;
        std::map<std::string, std::shared_ptr<Channel>> input_channels, output_channels;

        std::transform(streams.begin(), streams.end(), std::inserter(input_channels, input_channels.begin()),
                       [](auto &stream) { return std::make_pair(stream->key, std::make_shared<MessageChannel>()); }
        );

        std::transform(streams.begin(), streams.end(), std::inserter(output_channels, output_channels.begin()),
                       [](auto &stream) { return std::make_pair(stream->key, std::make_shared<MessageChannel>()); }
        );

        threads.emplace_back(error_handler.run(
                "parallel",
                [&]() { branch->process(input, input_channels, output, nested_handler); }
        ));

        threads.emplace_back(error_handler.run(
                "parallel",
                [&]() { merge->process(output_channels, output, nested_handler); }
        ));

        for (auto &stream : streams) {
            threads.emplace_back(error_handler.run(
                    "parallel",
                    [&](auto in, auto out) { stream->process(in, out, nested_handler); },
                    input_channels.at(stream->key),
                    output_channels.at(stream->key)
            ));
        }

        for(auto &thread : threads) { thread.join(); }
    }

    void Parallel::DecoratedBranch::process(
            std::shared_ptr<Channel> input,
            std::map<std::string, std::shared_ptr<Channel>> output,
            std::shared_ptr<Channel> bypass,
            ErrorHandler &error_handler
    ) {
        CloseGuard input_closer{input}, output_closer{output};
        error_handler.handle(key, [&]() { branch->process(input, output, bypass); });
    }

    Parallel::DecoratedBranch::DecoratedBranch(
            std::unique_ptr<Branch> branch,
            std::string key
    ) : branch(std::move(branch)), key(std::move(key)) {}

    void Parallel::DecoratedMerge::process(
            std::map<std::string, std::shared_ptr<Channel>> input,
            std::shared_ptr<Channel> output,
            ErrorHandler &error_handler
    ) {
        CloseGuard input_closer{input}, output_closer{output};
        error_handler.handle(key, [&]() { merge->process(input, output); });
    }

    Parallel::DecoratedMerge::DecoratedMerge(
            std::unique_ptr<Merge> merge,
            std::string key
    ) : merge(std::move(merge)), key(std::move(key)) {}
}