//
// Created by dchansen on 2/18/19.
//

#include "PureDistributed.h"
#include "connection/distributed/RemoteChannel.h"
#include "connection/distributed/remote_workers.h"

namespace {
    using namespace Gadgetron::Server::Distributed;
    using namespace Gadgetron::Core;

    class RemoteWorker {
    public:
        void process(InputChannel& input, OutputChannel& output) {
            for (auto message : input){
                try {
                    remote.push_message(message.clone());
                    output.push_message(remote.pop());
                } catch(const RemoteError& error) {
                    GDEBUG_STREAM(error.what())
                    retry_failed(std::move(message));
                    return;
                }
            }
        }

    private:
        RemoteChannel remote;
        std::function<void(Message)> retry_failed;
    };

}
void Gadgetron::Server::Connection::Stream::PureDistributed::process(Gadgetron::Core::InputChannel input,
    Gadgetron::Core::OutputChannel output, Gadgetron::Server::Connection::ErrorHandler& error_handler) {


}

const std::string& Gadgetron::Server::Connection::Stream::PureDistributed::name() {
    const static std::string n = "PureDistributed";
    return n;
}

Gadgetron::Server::Connection::Stream::PureDistributed::PureDistributed(
    const Gadgetron::Server::Connection::Config::PureDistributed& config, const Gadgetron::Core::Context& context, Loader& loader): loader{loader}, context{context} {
    readers = loader.load_readers(config);
    writers = loader.load_writers(config);
    remote_config = Connection::Config{config.readers,config.writers,{"PureStream",config.stream.gadgets}};
}

