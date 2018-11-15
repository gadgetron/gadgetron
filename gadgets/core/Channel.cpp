#include "Channel.h"

namespace Gadgetron {
    namespace Core {

        std::unique_ptr<Message>&& MessageChannel::pop() {
            std::unique_lock<std::mutex> lock(m);
            cv.wait(lock,[this](){return !this->queue.empty();});
            auto ptr = std::move(queue.front());
            queue.pop_front();
            return std::move(ptr);

        }


        void MessageChannel::push_message(std::unique_ptr<Message> && message ) {
            std::unique_lock<std::mutex> lock(m);
            queue.emplace_back(message);
            cv.notify_one();
        }


    }
}