#pragma once

#include <chrono>
#include <memory>
#include <future>

#include "connection/Core.h"
#include "connection/stream/common/Serialization.h"
#include "connection/stream/common/Configuration.h"
#include "connection/stream/common/ExternalChannel.h"

#include "Discovery.h"
#include "Pool.h"

#include "Message.h"

namespace Gadgetron::Server::Connection::Stream {

    class Worker {
    public:
        const Address address;

        Worker(
                const Address &address,
                std::shared_ptr<Serialization> serialization,
                std::shared_ptr<Configuration> configuration
        );

        void on_load_change(std::function<void()> callback);
        long long current_load() const;

        std::future<Core::Message> push(Core::Message message);

        std::thread start(ErrorHandler &);
        void close();

    private:
        mutable std::mutex mutex;
        std::unique_ptr<ExternalChannel> channel;

        struct Job {
            std::chrono::time_point<std::chrono::steady_clock> start;
            std::promise<Core::Message> response;
        };

        std::list<Job> jobs;
        std::list<std::function<void()>> callbacks;

        struct {
            std::chrono::milliseconds latest = std::chrono::seconds(5);
        } timing;

        void load_changed();
        void handle_inbound_messages();
        void process_inbound_message(Core::Message message);
    };

    std::unique_ptr<Worker> create_worker(
            Address address,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration
    );
}

