#include "Channel.h"

namespace Gadgetron {
    namespace Core {

        Message MessageChannel::pop() {
            std::unique_lock<std::mutex> lock(m);
            cv.wait(lock,[this](){return !this->queue.empty() || closed;});
            if (this->queue.empty()){
                throw ChannelClosed();
            }
            Message message = std::move(queue.front());
            queue.pop_front();
            return message ;

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
}