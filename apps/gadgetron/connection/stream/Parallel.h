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
        using Branch = Core::Parallel::Branch;
        using Merge  = Core::Parallel::Merge;

    public:
        Parallel(const Config::Parallel &, const Core::Context &, Loader &);

        void process(
                std::shared_ptr<Core::Channel> input,
                std::shared_ptr<Core::Channel> output,
                ErrorHandler &error_handler
        ) override;


        class DecoratedBranch {
        public:
            DecoratedBranch(std::unique_ptr<Branch> branch, std::string key);

            void process(
                    std::shared_ptr<Channel> input,
                    std::map<std::string, std::shared_ptr<Channel>> output,
                    std::shared_ptr<Channel> bypass,
                    ErrorHandler &error_handler
            );
        private:
            const std::string key;
            std::unique_ptr<Branch> branch;
        };

        class DecoratedMerge {
        public:
            DecoratedMerge(std::unique_ptr<Merge> merge, std::string key);

            void process(
                    std::map<std::string, std::shared_ptr<Channel>> input,
                    std::shared_ptr<Channel> output,
                    ErrorHandler &error_handler
            );
        private:
            const std::string key;
            std::unique_ptr<Merge> merge;
        };

    private:
        std::unique_ptr<DecoratedBranch> branch;
        std::unique_ptr<DecoratedMerge>  merge;
        std::vector<std::unique_ptr<Stream>> streams;
    };
}