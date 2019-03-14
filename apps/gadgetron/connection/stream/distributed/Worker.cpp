
#include "Worker.h"

#include <chrono>
#include <future>

#include "connection/stream/common/External.h"
#include "connection/stream/common/ExternalChannel.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection::Stream;

namespace {

    Remote as_remote(Local, const std::shared_ptr<Configuration> &configuration) {
        return Remote {
            "localhost",
            std::to_string(configuration->context.args["port"].as<unsigned short>())
        };
    }

    Remote as_remote(Remote remote, const std::shared_ptr<Configuration> &) {
        return remote;
    }

    template<class T>
    std::chrono::milliseconds time_since(std::chrono::time_point<T> instance) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() -
                instance
        );
    }

}

namespace Gadgetron::Server::Connection::Stream {

    Worker::Worker(
            const Address &address,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration
    ) : address(address) {
        channel = std::make_unique<ExternalChannel>(
                connect(boost::apply_visitor(
                        [&](auto address) { return as_remote(address, configuration); },
                        address
                )),
                std::move(serialization),
                std::move(configuration)
        );
    }

    void Worker::load_changed() { for (auto &f : callbacks) { f(); } }

    void Worker::on_load_change(std::function<void()> callback) {
        callbacks.push_back(std::move(callback));
    }

    long long Worker::current_load() const {
        std::lock_guard guard(mutex);

        auto current_job_duration_estimate = std::max(
                timing.latest.count(),
                time_since(jobs.front().start).count()
        );

        return current_job_duration_estimate * jobs.size();
    }

    std::future<Message> Worker::push(Message message) {

        Job job {
            std::chrono::steady_clock::now(),
            std::promise<Message>()
        };

        auto future = job.response.get_future();

        {
            std::lock_guard guard(mutex);
            channel->push_message(std::move(message));
            jobs.push_back(std::move(job));
        }

        load_changed();

        return future;
    }


    std::thread Worker::start(ErrorHandler &error_handler) {

        ErrorHandler nested_handler{
            error_handler,
            boost::apply_visitor([](auto a) { return to_string(a); }, address)
        };

        return nested_handler.run([=]() { handle_inbound_messages(); });
    }

    void Worker::close() {
        channel->close();
    }

    void Worker::handle_inbound_messages() {

        while(true) {
            try {
                process_inbound_message(channel->pop());
            }
            catch (ChannelClosed &) {
                break;
            }
            catch (std::exception &e) {
                // Things!
            }
        }
    }

    void Worker::process_inbound_message(Core::Message message) {

        Job job;
        {
            std::lock_guard guard(mutex);

            job = std::move(jobs.front()); jobs.pop_front();
            timing.latest = time_since(job.start);
        }

        load_changed();

        job.response.set_value(std::move(message));
    }
}