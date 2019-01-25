#include "Channel.h"


namespace Gadgetron::Core {

    Message MessageChannel::pop_impl(std::unique_lock<std::mutex> lock) {
        cv.wait(lock,[this](){return !this->queue.empty() || closed;});
        if (this->queue.empty()){
            throw ChannelClosed();
        }
        Message message = std::move(queue.front());
        queue.pop_front();
        return message ;
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

    void MessageChannel::push_message(Message message ) {
        std::unique_lock<std::mutex> lock(m);
        if (closed){
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
}
