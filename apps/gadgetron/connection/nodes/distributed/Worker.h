#pragma once

#include <chrono>
#include <memory>
#include <future>

#include "connection/Core.h"
#include "connection/nodes/common/Serialization.h"
#include "connection/nodes/common/Configuration.h"
#include "connection/nodes/common/ExternalChannel.h"
#include "connection/nodes/common/Discovery.h"

#include "Message.h"

namespace Gadgetron::Server::Connection::Nodes {

    class Worker {
    public:
        const Address address;

        ~Worker();
        Worker(
                Address address,
                std::shared_ptr<Serialization> serialization,
                std::shared_ptr<Configuration> configuration
        );

        std::future<Core::Message> push(Core::Message message);
        long long current_load() const;
        void close();

    private:
        mutable std::mutex mutex;

        std::thread inbound_thread;

        struct Timing {
            std::chrono::milliseconds latest = std::chrono::seconds(5);
        } timing;

        struct Job;
        std::list<Job> jobs;
        std::unique_ptr<ExternalChannel> channel;

        struct PushModule; struct LoadModule; struct ClosedPushModule; struct ClosedLoadModule;
        std::unique_ptr<PushModule> push_module;
        std::unique_ptr<LoadModule> load_module;
        void switch_to_closed_modules();

        void handle_inbound_messages();
        void process_inbound_message(Core::Message message);
        void fail_pending_messages(const std::exception_ptr &e);
    };
}
