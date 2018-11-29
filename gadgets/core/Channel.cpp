#include "Channel.h"

namespace Gadgetron {
    namespace Core {

        std::unique_ptr<Message> MessageChannel::pop() {
            std::unique_lock<std::mutex> lock(m);
            cv.wait(lock,[this](){return !this->queue.empty() || closed;});
            if (this->queue.empty()){
                throw ChannelClosedError();
            }
            std::unique_ptr<Message> ptr = std::move(queue.front());
            queue.pop_front();
            return ptr;

        }


        void MessageChannel::push_message(std::unique_ptr<Message> && message ) {
            std::unique_lock<std::mutex> lock(m);
            if (closed){
                throw ChannelClosedError();
            }
            queue.emplace_back(std::move(message));
            cv.notify_one();
        }

        void MessageChannel::close() {
            std::unique_lock<std::mutex> lock(m);
            closed = true;
            cv.notify_all();
        }


        void OutputChannel::push(std::unique_ptr<Gadgetron::Core::Message> && message) {
            this->push_message(std::move(message));
        }


    }
}