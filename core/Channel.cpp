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

    Message GenericInputChannel::pop() {
        return channel->pop();
    }

    optional<Message> GenericInputChannel::try_pop() {
        return channel->try_pop();
    }

    GenericInputChannel::GenericInputChannel(std::shared_ptr<Channel> channel) : channel{channel},
                                                                   closer{std::make_shared<Channel::Closer>(channel)} {
    }

    void OutputChannel::push_message(Gadgetron::Core::Message message) {
        channel->push_message(std::move(message));
    }

    OutputChannel::OutputChannel(std::shared_ptr<Channel> channel) : channel{channel},
                                                                     closer{std::make_shared<Channel::Closer>(
                                                                             channel)} {}

    ChannelIterator<OutputChannel> OutputChannel::begin() {
        return ChannelIterator<OutputChannel>(this);
    }
}


Gadgetron::Core::GenericInputChannel Gadgetron::Core::split(const GenericInputChannel&channel) {
    return GenericInputChannel{channel};
}

Gadgetron::Core::OutputChannel Gadgetron::Core::split(const OutputChannel &channel) {
    return OutputChannel{channel};
}
