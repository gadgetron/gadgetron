
#include "PureDistributed.h"

#include "distributed/Discovery.h"
#include "distributed/Worker.h"
#include "distributed/Pool.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection::Stream;

namespace {

    void add_worker_to_pool(
            Address address,
            const std::shared_ptr<Serialization> serialization,
            const std::shared_ptr<Configuration> configuration,
            std::shared_ptr<Pool<Worker>> workers
    ) {
        try {
            workers->add(std::make_unique<Worker>(address, serialization, configuration));
        }
        catch (std::exception &e) {
            GWARN_STREAM("Failed to initialize worker with address: " << address << " (" << e.what() << ")");
        }
    }
}

namespace Gadgetron::Server::Connection::Stream {

    class PureDistributed::Job {
    public:
        Message message;
        std::shared_ptr<Worker> worker;
        std::future<Message> response;
    };

    PureDistributed::Job PureDistributed::send_message_to_worker(Message message, std::shared_ptr<Worker> worker) {
        return Job {
            message.clone(),
            worker,
            worker->push(std::move(message))
        };
    }

    Message PureDistributed::get_message_from_worker(Job job, size_t retries) {
        try {
            return job.response.get();
        }
        catch (std::exception &e) {
            GWARN_STREAM("Worker " << job.worker->address << " reported error: " << e.what());

            if (!retries) throw std::runtime_error("Multiple workers failed processing job. Assuming terminal problem. Closing stream.");

            auto worker = workers->best();
            GWARN_STREAM("Job will be retried on worker: " << worker->address << " (" << retries << " retries left)");

            return get_message_from_worker(
                    send_message_to_worker(
                        std::move(job.message),
                        worker
                    ),
                    retries - 1
            );
        }
    }

    void PureDistributed::process_outbound(InputChannel input, Queue &jobs) {
        for (auto message : input) {
            jobs.push(send_message_to_worker(
                    std::move(message),
                    workers->best()
            ));
        }
        workers->close();
        jobs.close();
    }

    void PureDistributed::process_inbound(OutputChannel output, Queue &jobs) {
        while(true) {
            output.push(get_message_from_worker(jobs.pop()));
        }
    }

    void PureDistributed::process(
            InputChannel input,
            OutputChannel output,
            ErrorHandler& error_handler
    ) {
        auto queue = Queue();

        auto outbound = error_handler.run(
                [&](auto input) { process_outbound(std::move(input), queue); },
                std::move(input)
        );

        auto inbound = error_handler.run(
                [&](auto output) { process_inbound(std::move(output), queue); },
                std::move(output)
        );

        outbound.join(); inbound.join();
    }

    PureDistributed::PureDistributed(
            const Config::PureDistributed& config,
            const Context& context,
            Loader& loader
    ) : serialization(std::make_shared<Serialization>(
                loader.load_readers(config),
                loader.load_writers(config)
        )),
        configuration(std::make_shared<Configuration>(
                context,
                config
        )),
        workers{
            std::make_shared<Pool<Worker>>()
        }{
        trigger_node_discovery();
    }

    void PureDistributed::trigger_node_discovery() {
        for (const auto &address : discover_peers()) {
            std::thread{
                [=](auto... args) { add_worker_to_pool(args...); },
                address,
                serialization,
                configuration,
                workers
            }.detach();
        }
    }

    const std::string& Gadgetron::Server::Connection::Stream::PureDistributed::name() {
        const static std::string n = "PureDistributed";
        return n;
    }
}
