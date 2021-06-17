
#include "Pool.h"

#include "connection/nodes/common/Discovery.h"

#include "log.h"
#include "io/iostream_operators.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Server::Connection::Nodes;
using namespace Gadgetron;

namespace {

    bool compare_load(const std::unique_ptr<Worker> &a, const std::unique_ptr<Worker> &b) {
        return a->current_load() < b->current_load();
    }

    std::unique_ptr<Worker> &select_best_worker(std::list<std::unique_ptr<Worker>> &workers) {
        return *std::min_element(workers.begin(), workers.end(), compare_load);
    }

    Message send_to_best_worker(
            Message message,
            std::shared_ptr<std::list<std::unique_ptr<Worker>>> workers,
            size_t retries
    ) {

        if (!retries) throw std::runtime_error("Multiple workers failed processing job; aborting.");

        auto &worker = select_best_worker(*workers);

        try {
            auto response = worker->push(message.clone());
            GDEBUG_STREAM("Pushed message; waiting for response from worker " << worker->address);
            auto response_message = response.get();
            GDEBUG_STREAM("Response gotten from worker " << worker->address);
            return std::move(response_message);
        }
        catch (const std::exception &e) {
            GWARN_STREAM("Worker " << worker->address << " failed processing job. The job will be retried. [" << e.what() << "]");
            return send_to_best_worker(std::move(message), std::move(workers), retries - 1);
        }
    }
}

namespace Gadgetron::Server::Connection::Nodes {

    Pool::Pool(
            std::list<std::unique_ptr<Worker>> workers
    ) : workers(std::make_shared<std::list<std::unique_ptr<Worker>>>(std::move(workers))) {}

    std::future<Message> Pool::push(Message message) {
        return std::async(
                std::launch::async,
                send_to_best_worker,
                std::move(message),
                workers,
                3
        );
    }
}

