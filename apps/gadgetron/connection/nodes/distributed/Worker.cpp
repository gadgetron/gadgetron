
#include "Worker.h"

#include <chrono>
#include <future>

#include "connection/nodes/common/External.h"
#include "connection/nodes/common/ExternalChannel.h"

#include "io/iostream_operators.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection::Nodes;

namespace {

    template<class T>
    std::chrono::milliseconds time_since(std::chrono::time_point<T> instance) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - instance
        );
    }
}

namespace Gadgetron::Server::Connection::Nodes {

    struct Worker::Job {
        std::chrono::time_point<std::chrono::steady_clock> start;
        std::promise<Core::Message> response;
    };

    struct Module {
        Worker& worker;
        explicit Module(Worker &worker) : worker(worker) {}
        virtual ~Module() = default;
    };

    struct Worker::PushModule : public Module {
        using Module::Module;

        virtual std::future<Message> push(Message message) {
            GDEBUG_STREAM("Pushing message to remote worker " << worker.address);

            Job job {
                    std::chrono::steady_clock::now(),
                    std::promise<Message>()
            };

            auto future = job.response.get_future();
            worker.channel->push_message(std::move(message));
            worker.jobs.push_back(std::move(job));
            return future;
        };
    };

    struct Worker::LoadModule : public Module {
        using Module::Module;

        virtual long long current_load() {
            auto current_job_duration_estimate = std::max(
                    worker.timing.latest.count(),
                    time_since(worker.jobs.front().start).count()
            );

            return current_job_duration_estimate * worker.jobs.size();
        };
    };

    struct Worker::ClosedPushModule : public Worker::PushModule {
        using Worker::PushModule::PushModule;
        std::future<Message> push(Message message) override {
            throw std::runtime_error("Cannot push message to closed/failed worker.");
        }
    };

    struct Worker::ClosedLoadModule : public Worker::LoadModule {
        using Worker::LoadModule::LoadModule;
        long long current_load() override { return std::numeric_limits<long long>::max(); }
    };
}


namespace Gadgetron::Server::Connection::Nodes {

    Worker::~Worker() {
        channel->close();
        inbound_thread.join();
    }

    Worker::Worker(
            Address address,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration
    ) : address(std::move(address)) {
        GDEBUG_STREAM("Creating worker " << this->address);

        channel = std::make_unique<ExternalChannel>(
                connect(this->address, configuration),
                std::move(serialization),
                std::move(configuration)
        );

        load_module = std::make_unique<LoadModule>(*this);
        push_module = std::make_unique<PushModule>(*this);

        inbound_thread = std::thread([=]() { handle_inbound_messages(); });
    }

    long long Worker::current_load() const {
        std::lock_guard<std::mutex> guard(mutex);
        return load_module->current_load();
    }

    std::future<Message> Worker::push(Message message) {
        std::lock_guard<std::mutex> guard(mutex);
        return push_module->push(std::move(message));
    }

    void Worker::close() {
        std::lock_guard<std::mutex> guard(mutex);
        channel->close();
    }

    void Worker::handle_inbound_messages() {
        try {
            while(true) process_inbound_message(channel->pop());
        }
        catch (const ChannelClosed &) {}
        catch (const std::exception &e) {
            GWARN_STREAM("Worker " << address << " failed: " << e.what());
            fail_pending_messages(std::current_exception());
        }
        switch_to_closed_modules();
    }

    void Worker::process_inbound_message(Core::Message message) {
        GDEBUG_STREAM("Received message from remote worker " << address);

        std::lock_guard<std::mutex> guard(mutex);

        auto job = std::move(jobs.front()); jobs.pop_front();
        timing.latest = time_since(job.start);
        job.response.set_value(std::move(message));
    }

    void Worker::fail_pending_messages(const std::exception_ptr &e) {
        std::lock_guard<std::mutex> guard(mutex);

        for (auto &job : jobs) {
            job.response.set_exception(e);
        }
    }

    void Worker::switch_to_closed_modules() {
        std::lock_guard<std::mutex> guard(mutex);
        load_module = std::make_unique<ClosedLoadModule>(*this);
        push_module = std::make_unique<ClosedPushModule>(*this);
    }
}
