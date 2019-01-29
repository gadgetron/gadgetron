//
// Created by dchansen on 1/16/19.
//

#include "Distributed.h"
#include <algorithm>

namespace {
    using Worker = Gadgetron::Server::Connection::Stream::Distributed::Worker;
    using Address = Gadgetron::Server::Distributed::Address;
    using Local = Gadgetron::Server::Distributed::Local;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron;

    std::vector<Worker> get_workers() {
        return {Address{"localhost", "9003"},
                Address{"localhost", "9004"},
                Local{}};
    }


    class WorkerChannelCreator : public Gadgetron::Core::Distributed::ChannelCreator {
    public:
        WorkerChannelCreator(const Config::Distributed &distributed_config,
                             Gadgetron::Core::Context context,
                             Gadgetron::Server::Connection::Loader &loader,
                             std::shared_ptr<Core::Channel> output_channel,
                             ErrorHandler &error_handler)
                : loader(loader), output_channel{std::move(output_channel)}, context{std::move(context)},
                  xml_config{serialize_config(
                          Config{distributed_config.readers,
                                 distributed_config.writers,
                                 distributed_config.stream})},
                  workers{get_workers()}, error_handler{error_handler},
                  current_worker{make_cyclic(workers.begin(),
                                             workers.end())} {

            readers = loader.load_readers(distributed_config);
            writers = loader.load_writers(distributed_config);
        }

        std::shared_ptr<Core::OutputChannel> create() override {

            auto previous_worker = current_worker;
            while (true) {
                try {

                    auto result = boost::apply_visitor([&](auto v) { return this->create_channel(v); },
                                                       *current_worker);

                    ++current_worker;
                    return result;
                } catch (const std::runtime_error &) {
                    ++current_worker;
                    if (current_worker == previous_worker) throw;
                }
            }
        }

        ~WorkerChannelCreator() override {
            for (auto &c : channels)
                c->close();
            for (auto &wt : worker_threads)
                wt.join();


        }


    private:
        std::shared_ptr<Core::Channel> create_channel(Local) {
            auto result = std::make_shared<Core::MessageChannel>();

            worker_threads.emplace_back(error_handler.run("Distributed local stream",
                                                          [&](auto input, auto output, auto config) {
                                                              auto stream = Stream::Stream{config, context, loader};
                                                              stream.process(input, output, error_handler);
                                                          }, result, output_channel, stream_config));
            channels.push_back(result);
            return result;
        }


        std::shared_ptr<Core::Channel> create_channel(Address address) {
            auto result = std::make_shared<Gadgetron::Server::Distributed::RemoteChannel>(address, xml_config,
                                                                                          context.header, readers,
                                                                                          writers);

            worker_threads.emplace_back(error_handler.run("RemoteChannel reader",
                                                          [](auto input, auto output) {

                                                              Core::InputChannel &inputChannel = *input;
                                                              Core::OutputChannel &outputChannel = *output;
                                                              std::transform(begin(inputChannel), end(inputChannel),
                                                                             begin(outputChannel),
                                                                             [](auto &&m) { return std::move(m); });
                                                          }, result, output_channel));
            channels.push_back(result);
            return result;
        }

        const std::string xml_config;
        const Config::Stream stream_config;
        const Core::Context context;
        const std::vector<Worker> workers;

        CyclicIterator<decltype(workers)::const_iterator> current_worker;
        std::map<uint16_t, std::unique_ptr<Core::Reader>> readers;
        std::vector<std::unique_ptr<Core::Writer>> writers;
        Gadgetron::Server::Connection::Loader &loader;

        std::vector<std::thread> worker_threads;
        std::vector<std::shared_ptr<Core::Channel>> channels;
        std::shared_ptr<Core::Channel> output_channel;
        ErrorHandler &error_handler;
    };


    void
    process_channel(std::shared_ptr<Gadgetron::Core::Channel> input, std::shared_ptr<Gadgetron::Core::Channel> output) {

        Gadgetron::Core::InputChannel &in_view = *input;
        for (auto &&message : in_view)
            output->push_message(std::move(message));
    }


}

void Gadgetron::Server::Connection::Stream::Distributed::process(std::shared_ptr<Gadgetron::Core::Channel> input,
                                                                 std::shared_ptr<Gadgetron::Core::Channel> output,
                                                                 Gadgetron::Server::Connection::ErrorHandler &error_handler) {

    auto channelcreator = WorkerChannelCreator{config, context, loader, output, error_handler};

    error_handler.handle("Distributor", [&]() { distributor->process(*input, channelcreator, *output); });

    input->close();
    output->close();

}

Gadgetron::Server::Connection::Stream::Distributed::Distributed(const Config::Distributed &distributed_config,
                                                                const Gadgetron::Core::Context &context,
                                                                Gadgetron::Server::Connection::Loader &loader)
        : context{context}, loader{loader}, config(distributed_config) {
    distributor = load_distributor(distributed_config.distributor);


}

std::unique_ptr<Gadgetron::Core::Distributed::Distributor>
Gadgetron::Server::Connection::Stream::Distributed::load_distributor(
        const Gadgetron::Server::Connection::Config::Distributor &conf) {
    auto factory = loader.load_factory<Loader::generic_factory<Core::Distributed::Distributor>>(
            "distributor_factory_export_", conf.classname, conf.dll);
    return factory(context, conf.properties);
}


