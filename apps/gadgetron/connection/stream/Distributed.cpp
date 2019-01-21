//
// Created by dchansen on 1/16/19.
//

#include "Distributed.h"


void Gadgetron::Server::Connection::Stream::Distributed::process(std::shared_ptr<Gadgetron::Core::Channel> input,
                                                                 std::shared_ptr<Gadgetron::Core::Channel> output,
                                                                 Gadgetron::Server::Connection::ErrorHandler &error_handler) {

}

Gadgetron::Server::Connection::Stream::Distributed::Distributed(const Config::Distributed &distributed_config,
                                                                const Gadgetron::Core::Context &context,
                                                                Gadgetron::Server::Connection::Loader &loader)
        : stream_config(distributed_config.stream), loader(loader), context(context) {
    distributor  = load_distributor(distributed_config.distributor);

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
