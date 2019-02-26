//
// Created by dchansen on 2/18/19.
//

#include "PureDistributed.h"
#include "MPMCChannel.h"
#include "connection/distributed/RemoteChannel.h"
#include "connection/distributed/remote_workers.h"

namespace {
    using namespace Gadgetron::Server::Distributed;
    using namespace Gadgetron::Core;
    class RemoteWorker;
    using WorkerQueue = MPMCChannel<std::shared_ptr<RemoteWorker>>;

    class RemoteWorker : public std::enable_shared_from_this<RemoteWorker> {
    public:
        RemoteWorker(const Address& address, const std::string& xml_config, const ISMRMRD::IsmrmrdHeader& header,
            const std::map<uint16_t, std::unique_ptr<Gadgetron::Core::Reader>>& readers,
            const std::vector<std::unique_ptr<Gadgetron::Core::Writer>>& writers, WorkerQueue& worker_queue)
            : worker_queue{ worker_queue } {
            auto channels = make_channel<RemoteChannel>(address, xml_config, header, readers, writers);

            input_writer
                = std::thread([this](OutputChannel out) { write_input(std::move(out)); }, std::move(channels.output));

            output_reader
                = std::thread([this](InputChannel in) { read_output((std::move(in))); }, std::move(channels.input));
        }

        std::future<Message> remote_process(Message message) {
            std::lock_guard guard(write_lock);
            work.push(std::move(message));
            return get_output();
        }

        void close() {
            work.close();
            output_reader.join();
            input_writer.join();
        }

    protected:
        std::future<Message> get_output() {
            outputs.emplace_back();
            return outputs.back().get_future();
        }

        void read_output(InputChannel remoteInput) {
            try {
                for (auto message : remoteInput) {

                    std::lock_guard guard(write_lock);
                    outputs.front().set_value(std::move(message));
                    outputs.pop_front();
                    worker_queue.push(this->shared_from_this());
                }
            } catch (...) {
                handle_error(std::current_exception());
            }
        }

        void write_input(OutputChannel remoteOutput) {
            try {
                for (;;) {
                    remoteOutput.push_message(work.pop());
                }
            } catch (const ChannelClosed&) {
            } catch (...) {
                handle_error(std::current_exception());
            }
        }

        void handle_error(std::exception_ptr exception) {
            std::lock_guard guard(write_lock);
            for (auto& promises : outputs) {
                promises.set_exception(exception);
            }
            outputs.clear();
        }

    private:
        WorkerQueue& worker_queue;
        std::mutex write_lock;
        std::list<std::promise<Message>> outputs;
        MPMCChannel<Message> work;

        std::thread output_reader;
        std::thread input_writer;
    };

    template <class Container> Container take_n(const Container& container, size_t n) {
        Container c;
        c.reserve(n);
        for (size_t i = 0; i < n && i < container.size(); i++)
            c.push_back(container[i]);
        return c;
    }



    //
    class OutputHandler {
    public:
        OutputHandler(WorkerQueue& worker_queue, OutputChannel output, Gadgetron::Server::Connection::ErrorHandler& error_handler)
            : worker_queue{ worker_queue } {
            writer = error_handler.run([this](auto channel) { this->check_futures(std::move(channel)); },std::move(output));
        }

        void add_work(std::future<Message> future_message, Message message) {
            work.emplace(std::move(future_message), std::move(message));
        }

        void close(){
            work.close();
            writer.join();
        }

    private:
        void check_futures(OutputChannel output) {

            for (;;) {
                auto done_work = work.pop();
                try {
                    auto message = done_work.first.get();
                    output.push_message(std::move(message));
                } catch (const RemoteError& remoteerror) {
                    throw;
                } catch (...){
                    auto future_work = worker_queue.pop()->remote_process(done_work.second.clone());
                    work.emplace(std::move(future_work),std::move(done_work.second));
                }

            }
        }

        std::thread writer;
        WorkerQueue& worker_queue;
        MPMCChannel<std::pair<std::future<Message>, Message>> work;
    };

}

void Gadgetron::Server::Connection::Stream::PureDistributed::process(Gadgetron::Core::InputChannel input,
    Gadgetron::Core::OutputChannel output, Gadgetron::Server::Connection::ErrorHandler& error_handler) {

    auto available_workers = WorkerQueue();
    auto remotes           = nworkers ? take_n(get_remote_workers(), nworkers) : get_remote_workers();
    auto workers           = std::vector<std::shared_ptr<RemoteWorker>>();
    std::transform(remotes.begin(), remotes.end(), std::back_inserter(workers), [&](const auto& address) {
        return std::make_shared<RemoteWorker>(
            address, remote_config, context.header, readers, writers, available_workers);
    });

    if (workers.empty()) throw std::runtime_error("No workers");
    for (int i = 0; i < 2; i++)
        for (auto& worker : workers)
            available_workers.push(worker);

    auto outputhandler = OutputHandler(available_workers,std::move(output),error_handler);

    for (auto message : input) {
        auto result = available_workers.pop()->remote_process(message.clone());
        outputhandler.add_work(std::move(result),std::move(message));
    }

    for (auto& worker : workers)
        worker->close();

    outputhandler.close();

}

const std::string& Gadgetron::Server::Connection::Stream::PureDistributed::name() {
    const static std::string n = "PureDistributed";
    return n;
}

Gadgetron::Server::Connection::Stream::PureDistributed::PureDistributed(
    const Gadgetron::Server::Connection::Config::PureDistributed& config, const Gadgetron::Core::Context& context,
    Loader& loader)
    : loader{ loader }, context{ context } {
    readers       = loader.load_readers(config);
    writers       = loader.load_writers(config);

    std::vector<Config::Node> nodes(config.stream.gadgets.begin(),config.stream.gadgets.end());
    remote_config = serialize_config(Connection::Config{ config.readers, config.writers,
        Config::Stream{
            "PureStream", nodes } });
}
