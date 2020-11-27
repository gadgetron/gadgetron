#pragma once

#include "connection/stream/Processable.h"
#include "connection/stream/Stream.h"

#include "parallel/Branch.h"
#include "parallel/Merge.h"

#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::Stream {

    class Parallel : public Processable {
        using Channel = Core::Channel;
        using InputChannel = Core::GenericInputChannel;
        using OutputChannel = Core::OutputChannel;
        using Branch = Core::Parallel::Branch;
        using Merge  = Core::Parallel::Merge;

    public:
        Parallel(const Config::Parallel &, const Core::StreamContext &, Loader &);

        void process(
                Core::GenericInputChannel input,
                Core::OutputChannel output,
                ErrorHandler &error_handler
        ) override;

        const std::string &name() override;

        class DecoratedBranch {
        public:
            DecoratedBranch(std::unique_ptr<Branch> branch, std::string key);

            const std::string key;
            void process(
                    InputChannel input,
                    std::map<std::string, OutputChannel> output,
                    OutputChannel bypass
            );
        private:
            std::unique_ptr<Branch> branch;
        };

        class DecoratedMerge {
        public:
            DecoratedMerge(std::unique_ptr<Merge> merge, std::string key);

            const std::string key;
            void process(
                    std::map<std::string, InputChannel> input,
                    OutputChannel output
            );
        private:
            std::unique_ptr<Merge> merge;
        };

    private:
        std::unique_ptr<DecoratedBranch> branch;
        std::unique_ptr<DecoratedMerge>  merge;
        std::vector<std::shared_ptr<Stream>> streams;
    };

}