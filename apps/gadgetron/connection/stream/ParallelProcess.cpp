//
// Created by dchansen on 2/7/19.
//

#include "ParallelProcess.h"
#include "omp.h"

void Gadgetron::Server::Connection::Stream::ParallelProcess::process(Gadgetron::Core::InputChannel input,
    Gadgetron::Core::OutputChannel output, Gadgetron::Server::Connection::ErrorHandler& error_handler) {

    const size_t n_threads = workers ? workers : omp_get_num_threads();
    GDEBUG("Starting parallel process with %d workers\n", n_threads);
    std::vector<std::thread> threads;
    for (size_t worker = 0; worker < n_threads; worker++) {
        threads.push_back(error_handler.run([&]() {
            try {
                for (;;) {
                    output.push_message(pureStream.process_function(input.pop()));
                }

            } catch (const Core::ChannelClosed&) {
            }
        }));
    }

    for (auto& t : threads)
        t.join();
}

Gadgetron::Server::Connection::Stream::ParallelProcess::ParallelProcess(
    const Gadgetron::Server::Connection::Config::ParallelProcess& conf, const Gadgetron::Core::Context& context,
    Gadgetron::Server::Connection::Loader& loader)
    : pureStream{ conf.stream, context, loader }, workers{ conf.workers } {}

const std::string& Gadgetron::Server::Connection::Stream::ParallelProcess::name() {
    const static std::string n = "ParallelProcess";
    return n;
}
