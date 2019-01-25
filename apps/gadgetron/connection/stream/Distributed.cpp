//
// Created by dchansen on 1/16/19.
//

#include "Distributed.h"

namespace {
    std::vector<Gadgetron::Server::Distributed::Address> get_workers() {
        return {{"localhost", "9003"}};
    }


    class ChannelCreatorFunction : public Gadgetron::Core::Distributed::ChannelCreator {
    public:
        template<class F>
        explicit ChannelCreatorFunction(F func) : func(func) {}

        std::shared_ptr<Gadgetron::Core::OutputChannel> create() override {
            return func();
        }

       std::function<std::shared_ptr<Gadgetron::Core::OutputChannel>(void)> func;

    };




    void
    process_channel(std::shared_ptr<Gadgetron::Core::Channel> input, std::shared_ptr<Gadgetron::Core::Channel> output) {

        Gadgetron::Core::InputChannel &in_view = *input;
        std::this_thread::sleep_for(std::chrono::seconds(10));
        for (auto &&message : in_view)
            output->push_message(std::move(message));
    }

    }

void Gadgetron::Server::Connection::Stream::Distributed::process(std::shared_ptr<Gadgetron::Core::Channel> input,
                                                                 std::shared_ptr<Gadgetron::Core::Channel> output,
                                                                 Gadgetron::Server::Connection::ErrorHandler &error_handler) {
    std::vector<std::thread> threads;
    std::vector<std::shared_ptr<RemoteChannel>> channels;

    auto creator = ChannelCreatorFunction([&]() {
        auto channel = this->create_remote_channel();

        threads.emplace_back(error_handler.run("ChannelReader", process_channel, channel, output));
        channels.push_back(channel);
        return channel;
    });

    error_handler.handle("Distributor",[&](){distributor->process(*input,creator,*output);});

    input->close();
    for (auto channel : channels)
        channel->close();

    for (auto& t : threads)
        t.join();

    output->close();

}

Gadgetron::Server::Connection::Stream::Distributed::Distributed(const Config::Distributed &distributed_config,
                                                                const Gadgetron::Core::Context &context,
                                                                Gadgetron::Server::Connection::Loader &loader)
        : loader(loader), context(context),
          xml_config{serialize_config(
                  Config{distributed_config.readers, distributed_config.writers, distributed_config.stream})},
          workers{get_workers()} {

    distributor = load_distributor(distributed_config.distributor);

    readers = loader.load_readers(distributed_config);
    writers = loader.load_writers(distributed_config);


}

std::unique_ptr<Gadgetron::Core::Distributed::Distributor>
Gadgetron::Server::Connection::Stream::Distributed::load_distributor(
        const Gadgetron::Server::Connection::Config::Distributor &conf) {
    auto factory = loader.load_factory<Loader::generic_factory<Core::Distributed::Distributor>>(
            "distributor_factory_export_", conf.classname, conf.dll);
    return factory(context, conf.properties);
}

std::shared_ptr<Gadgetron::Server::Distributed::RemoteChannel>
Gadgetron::Server::Connection::Stream::Distributed::create_remote_channel() {

    return std::make_shared<RemoteChannel>(workers.front(), xml_config,context.header, readers, writers);
}
