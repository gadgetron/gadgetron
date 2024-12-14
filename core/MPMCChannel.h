#pragma once

#include <condition_variable>
#include <list>
#include <mutex>
#include <optional>

namespace Gadgetron::Core {

    template <class T> class MPMCChannel {
    public:
        MPMCChannel() = default;
        MPMCChannel(MPMCChannel&&) noexcept;
        void push(T);

        template <class... ARGS> void emplace(ARGS&&... args);

        T pop();
        std::optional<T> try_pop();

        void close();

    private:
        T pop_impl(std::unique_lock<std::mutex> lock);
        std::list<T> queue;
        bool is_closed = false;
        std::mutex m;
        std::condition_variable cv;
    };

    class ChannelClosed : public std::runtime_error {
    public:
        ChannelClosed() : std::runtime_error("Channel was closed"){};
    };

    /** Implementation **/

    template <class T> T MPMCChannel<T>::pop_impl(std::unique_lock<std::mutex> lock) {
        cv.wait(lock, [this]() { return !this->queue.empty() || is_closed; });
        if (queue.empty()) {
            throw ChannelClosed();
        }
        T message = std::move(queue.front());
        queue.pop_front();
        return message;
    }

    template <class T> T MPMCChannel<T>::pop() {
        std::unique_lock<std::mutex> lock(m);
        return pop_impl(std::move(lock));
    }

    template <class T> std::optional<T> MPMCChannel<T>::try_pop() {
        std::unique_lock<std::mutex> lock(m);
        if (queue.empty()) {
            return std::nullopt;
        }
        return pop_impl(std::move(lock));
    }

    template <class T> void MPMCChannel<T>::push(T message) {
        {
            std::lock_guard<std::mutex> lock(m);
            if (is_closed)
                throw ChannelClosed();
            queue.emplace_back(std::move(message));
        }
        cv.notify_one();
    }

    template <class T> void MPMCChannel<T>::close() {
        {
            std::lock_guard<std::mutex> lock(m);
            is_closed = true;
        }
        cv.notify_all();
    }

    template <class T> template <class... ARGS> void MPMCChannel<T>::emplace(ARGS&&... args) {
        {
            std::lock_guard<std::mutex> guard(m);
            queue.emplace_back(std::forward<ARGS>(args)...);
        }
        cv.notify_one();
    }
    template <class T> MPMCChannel<T>::MPMCChannel(MPMCChannel&& other) noexcept {
        std::lock_guard<std::mutex> guard(other.m);
        this->queue = std::move(other.queue);
        this->is_closed = other.is_closed;
        other.is_closed = true;
    }

}
