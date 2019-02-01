#include "Channel.h"


namespace Gadgetron::Core {

    class Channel::Closer {
    public:
        explicit Closer(std::shared_ptr<Channel> channel) : channel{channel} {}

        ~Closer() { channel->close(); }

    private:
        std::shared_ptr<Channel> channel;
    };

    Message MessageChannel::pop_impl(std::unique_lock<std::mutex> lock) {
        cv.wait(lock, [this]() { return !this->queue.empty() || closed; });
        if (queue.empty()) {
            throw ChannelClosed();
        }
        Message message = std::move(queue.front());
        queue.pop_front();
        return message;
    }

    Message MessageChannel::pop() {
        std::unique_lock<std::mutex> lock(m);
        return pop_impl(std::move(lock));
    }

    optional<Message> MessageChannel::try_pop() {
        std::unique_lock<std::mutex> lock(m);
        if (queue.empty()) {
            return none;
        }
        return pop_impl(std::move(lock));
    }

    void MessageChannel::push_message(Message message) {
        std::unique_lock<std::mutex> lock(m);
        if (closed) {
            throw ChannelClosed();
        }
        queue.emplace_back(std::move(message));
        cv.notify_one();
    }

    void MessageChannel::close() {
        std::unique_lock<std::mutex> lock(m);
        closed = true;
        cv.notify_all();
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
