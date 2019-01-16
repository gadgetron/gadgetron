#pragma once

#include <map>
#include <memory>

#include "Channel.h"

namespace Gadgetron::Server::Connection {

    class NodeHandler {
    public:
        virtual void process(std::shared_ptr<Core::Channel> in, std::shared_ptr<Core::Channel> out) = 0;
        virtual ~NodeHandler() = default;
    };

    class BranchHandler {
    public:
        virtual void process(
                std::shared_ptr<Core::Channel> input,
                const std::map<std::string, std::shared_ptr<Core::Channel>> &output,
                std::shared_ptr<Core::Channel> bypass
        ) = 0;
        virtual ~BranchHandler() = default;
    };

    class MergeHandler {
    public:
        virtual void process(
                std::map<std::string, std::shared_ptr<Core::Channel>> input,
                std::shared_ptr<Core::Channel> output
        ) = 0;
        virtual ~MergeHandler() = default;
    };
}
