#include "Channel.h"


namespace Gadgetron::Core {

    class Channel::Closer {
    public:
        explicit Closer(std::shared_ptr<Channel> channel) : channel{std::move(channel)} {}

        ~Closer() { channel->close(); }

    private:
        std::shared_ptr<Channel> channel;
    };

    Message MessageChannel::pop() {
        return channel.pop();
    }

    optional<Message> MessageChannel::try_pop() {
       return channel.try_pop();
    }

    void MessageChannel::push_message(Message message) {
        channel.push(std::move(message));
    }

    void MessageChannel::close() {
       channel.close();
    }

    Message InputChannel::pop() {
        return channel->pop();
    }

    optional<Message> InputChannel::try_pop() {
        return channel->try_pop();
    }

    InputChannel::InputChannel(std::shared_ptr<Channel> channel) : channel{channel},
                                                                   closer{std::make_shared<Channel::Closer>(channel)} {
    }

    void OutputChannel::push_message(Gadgetron::Core::Message message) {
        channel->push_message(std::move(message));
    }

    OutputChannel::OutputChannel(std::shared_ptr<Channel> channel) : channel{channel},
                                                                     closer{std::make_shared<Channel::Closer>(
                                                                             channel)} {}
}


Gadgetron::Core::InputChannel Gadgetron::Core::split(const InputChannel &channel) {
    return InputChannel{channel};
}

Gadgetron::Core::OutputChannel Gadgetron::Core::split(const OutputChannel &channel) {
    return OutputChannel{channel};
}
