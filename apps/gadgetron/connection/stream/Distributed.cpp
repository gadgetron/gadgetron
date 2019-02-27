#include <algorithm>

#include "Distributed.h"

#include "connection/distributed/remote_workers.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Distributed;

    std::vector<Worker> get_workers() {

        auto workers = get_remote_workers();

        if (workers.empty()) {
            GWARN_STREAM(
                    "Remote worker list empty; adding local worker. " <<
                    "This machine will perform reconstructions. " <<
                    "This is probably not what you intended."
            )
            workers.emplace_back(Local{});
        }

        return workers;
    }

    std::string print_worker(const Address& address) {
        return address.ip + ":" + address.port;
    }

    std::string print_worker(const Local& local) {
        return "Local";
    }

    class WorkerChannelCreator : public Core::Distributed::ChannelCreator {
    public:
        WorkerChannelCreator(
                const Config::Distributed& distributed_config,
                Gadgetron::Core::Context context,
                Gadgetron::Server::Connection::Loader& loader,
                Core::OutputChannel output_channel,
                ErrorHandler& error_handler
        ) : loader { loader },
            output_channel{ std::move(output_channel) },
            context{ std::move(context) },
            xml_config{ serialize_config(
                  Config{ distributed_config.readers, distributed_config.writers, distributed_config.stream })
            },
            workers{ get_workers() },
            error_handler{ error_handler },
            current_worker{ make_cyclic(workers.begin(), workers.end()) },
            stream_config{ distributed_config.stream } {

            readers = loader.load_readers(distributed_config);
            writers = loader.load_writers(distributed_config);
        }

        Core::OutputChannel create() override {

            auto previous_worker = current_worker;
            while (true) {
                try {
                    auto worker = *current_worker;
                    GDEBUG_STREAM(boost::apply_visitor([](auto w) { return print_worker(w); }, worker))
                    auto result = boost::apply_visitor([&](auto v) { return this->create_channel(v); }, worker);

                    ++current_worker;
                    return result;
                } catch (const std::runtime_error&) {

                    ++current_worker;
                    if (current_worker == previous_worker)
                        throw;
                }
            }
        }

        ~WorkerChannelCreator() override {
            for (auto& wt : worker_threads) wt.join();
        }

    private:

        Core::OutputChannel create_channel(Local) {

            auto channel = Core::make_channel<Core::MessageChannel>();

            worker_threads.emplace_back(
                    Stream::Processable::process_async(
                            std::make_shared<Stream::Stream>(stream_config, context, loader),
                            std::move(channel.input),
                            Core::split(output_channel),
                            ErrorHandler{ error_handler, "Distributed (Local Stream)" }
                    )
            );

            return std::move(channel.output);
        }

        Core::OutputChannel create_channel(Address address) {
            auto channel = Core::make_channel<Gadgetron::Server::Distributed::RemoteChannel>(
                address, xml_config, context.header, readers, writers);

            worker_threads.emplace_back(
                    ErrorHandler{ error_handler, "Distributed (Remote Stream)" }.run(
                            [](auto input, auto output) {
                                std::transform(begin(input), end(input), begin(output), [](auto&& m) { return std::forward<decltype(m)>(m); });
                            },
                            std::move(channel.input),
                            Core::split(output_channel)
                    )
            );

            return std::move(channel.output);
        }

        const std::string xml_config;
        const Config::Stream stream_config;
        const Core::Context context;
        const std::vector<Worker> workers;

        CyclicIterator<decltype(workers)::const_iterator> current_worker;
        std::map<uint16_t, std::unique_ptr<Core::Reader>> readers;
        std::vector<std::unique_ptr<Core::Writer>> writers;
        Gadgetron::Server::Connection::Loader& loader;

        std::vector<std::thread> worker_threads;
        Core::OutputChannel output_channel;
        ErrorHandler& error_handler;
    };

}

namespace Gadgetron::Server::Connection::Stream {

    void Distributed::process(
            Core::InputChannel input,
            Core::OutputChannel output,
            ErrorHandler& error_handler
    ) {
        auto channel_creator = WorkerChannelCreator{
            config,
            context,
            loader,
            Core::split(output),
            error_handler
        };


        distributor->process(std::move(input), channel_creator, std::move(output));
    }

    Distributed::Distributed(
            const Config::Distributed& distributed_config,
            const Core::Context& context,
            Loader& loader
    ) : context{ context },
        loader{ loader },
        config(distributed_config) {
        distributor = load_distributor(distributed_config.distributor);
    }

    std::unique_ptr<Core::Distributed::Distributor> Distributed::load_distributor(
            const Config::Distributor& conf
    ) {
        auto factory = loader.load_factory<Loader::generic_factory<Core::Distributed::Distributor>>(
                "distributor_factory_export_", conf.classname, conf.dll);
        return factory(context, conf.properties);
    }

    const std::string& Distributed::name() {
        static const std::string n = "Distributed";
        return n;
    }
}


