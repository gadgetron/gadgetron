#include "ParallelProcess.h"

#include "ThreadPool.h"

using namespace Gadgetron::Core;

namespace Gadgetron::Server::Connection::Nodes {

    void ParallelProcess::process_input(GenericInputChannel input, Queue &queue) {

        ThreadPool pool(workers ? workers : std::thread::hardware_concurrency());

        for (auto message : input) {
            queue.push(
                pool.async(
                        [&](auto message) { return pureStream.process_function(std::move(message)); },
                        std::move(message)
                )
            );
        }

        pool.join(); queue.close();
    }

    void ParallelProcess::process_output(OutputChannel output, Queue &queue) {
        while(true) output.push_message(queue.pop().get());
    }

    void ParallelProcess::process(GenericInputChannel input,
            OutputChannel output,
            ErrorHandler& error_handler
    ) {
        Queue queue;

        auto input_thread = error_handler.run(
                [&](auto input) { this->process_input(std::move(input), queue); },
                std::move(input)
        );

        auto output_thread = error_handler.run(
                [&](auto output) { this->process_output(std::move(output), queue); },
                std::move(output)
        );

        input_thread.join(); output_thread.join();
    }

    ParallelProcess::ParallelProcess(
            const Config::ParallelProcess& conf,
            const Context& context,
            Loader& loader
    ) : pureStream{ conf.stream, context, loader }, workers{ conf.workers } {}

    const std::string& ParallelProcess::name() {
        const static std::string n = "ParallelProcess";
        return n;
    }
}


