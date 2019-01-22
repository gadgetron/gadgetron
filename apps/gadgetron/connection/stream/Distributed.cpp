//
// Created by dchansen on 1/16/19.
//

#include "Distributed.h"


namespace {
    std::vector<Gadgetron::Server::Distributed::Address> get_workers() {
        return {{"localhost", "9002"}};
    }
}

void Gadgetron::Server::Connection::Stream::Distributed::process(std::shared_ptr<Gadgetron::Core::Channel> input,
                                                                 std::shared_ptr<Gadgetron::Core::Channel> output,
                                                                 Gadgetron::Server::Connection::ErrorHandler &error_handler) {

}

Gadgetron::Server::Connection::Stream::Distributed::Distributed(const Config::Distributed &distributed_config,
                                                                const Gadgetron::Core::Context &context,
                                                                Gadgetron::Server::Connection::Loader &loader)
        : loader(loader), context(context),
          config{distributed_config.readers, distributed_config.writers, distributed_config.stream},
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


    return std::shared_ptr<Gadgetron::Server::Distributed::RemoteChannel>();
}
