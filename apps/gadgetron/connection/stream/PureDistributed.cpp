
#include "PureDistributed.h"

#include "connection/stream/common/Discovery.h"
#include "connection/stream/common/Closer.h"

#include "distributed/Worker.h"
#include "distributed/Pool.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Stream;


namespace {

    std::unique_ptr<Worker> connect_to_peer(
            Address address,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration
    ) {
        return std::make_unique<Worker>(std::move(address), std::move(serialization), std::move(configuration));
    }

    std::list<std::future<std::unique_ptr<Worker>>> begin_connecting_to_peers(
            std::future<std::vector<Address>> addresses,
            std::shared_ptr<Serialization> serialization,
            std::shared_ptr<Configuration> configuration
    ) {
        std::list<std::future<std::unique_ptr<Worker>>> workers;
        for (auto address : addresses.get()) {
            workers.emplace_back(std::async(std::launch::async, connect_to_peer, address, serialization, configuration));
        }
        return std::move(workers);
    }

    std::list<std::unique_ptr<Worker>> finish_connecting_to_peers(
            std::list<std::future<std::unique_ptr<Worker>>> pending_workers
    ) {
        std::list<std::unique_ptr<Worker>> workers;
        for (auto &pending_worker : pending_workers) {
            auto status = pending_worker.wait_for(std::chrono::seconds(5));
            if (status == std::future_status::ready) {
                try {
                    workers.emplace_back(pending_worker.get());
                }
                catch (std::exception &e) {
                    GWARN_STREAM("Failed connecting to worker with error: " << e.what());
                }
            }
        }

        if (workers.empty()) throw std::runtime_error("Could not connect to ANY workers; cannot distribute work.");

        return std::move(workers);
    }
}

namespace Gadgetron::Server::Connection::Stream {


    void PureDistributed::process_outbound(GenericInputChannel input, std::shared_ptr<Queue> jobs) {

        auto closer = make_closer(jobs);

        auto workers = Pool(finish_connecting_to_peers(std::move(pending_workers)));

        for (auto message : input) {
            jobs->push(workers.push(std::move(message)));
        }
    }

    void PureDistributed::process_inbound(OutputChannel output, std::shared_ptr<Queue> jobs) {
        while (true) {
            output.push_message(jobs->pop().get());
        }
    }

    void PureDistributed::process(GenericInputChannel input,
            OutputChannel output,
            ErrorHandler& error_handler
    ) {
        auto queue = std::make_shared<Queue>();

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
            const Core::StreamContext& context,
            Loader& loader
    ) : serialization(std::make_shared<Serialization>(
                loader.load_readers(config),
                loader.load_writers(config)
        )),
        configuration(std::make_shared<Configuration>(
                context,
                config
        )) {
        pending_workers = begin_connecting_to_peers(std::async(discover_peers), serialization, configuration);
    }


    const std::string& Gadgetron::Server::Connection::Stream::PureDistributed::name() {
        const static std::string n = "PureDistributed";
        return n;
    }
}
