#pragma once

#include <mutex>
#include <vector>
#include <memory>
#include <algorithm>

namespace Gadgetron::Server::Connection::Stream {

    template<class T>
    class Pool : public std::enable_shared_from_this<Pool<T>> {
    public:
        void add(std::unique_ptr<T>);
        void close();
        std::shared_ptr<T> best();

    private:
        void on_failure(const std::shared_ptr<T>);
        void on_load_change();

        std::mutex mutex;
        std::vector<std::shared_ptr<T>> ts;
    };

    template<class T>
    void Pool<T>::add(std::unique_ptr<T> t) {
        std::lock_guard<std::mutex> guard(mutex);

        std::shared_ptr<T> worker = std::move(t);
        auto self = this->shared_from_this();

        worker->on_failure([=]() { self->on_failure(worker); });
        worker->on_load_change([=]() { self->on_load_change(); });

        ts.push_back(std::move(worker));
    }

    template<class T>
    void Pool<T>::close() {
        std::lock_guard<std::mutex> guard(mutex);
        for (auto &t : ts) t->close();
    }

    template<class T>
    std::shared_ptr<T> Pool<T>::best() {
        std::lock_guard<std::mutex> guard(mutex);
        if (ts.empty()) throw std::runtime_error("Worker pool empty; no best worker available.");
        return *ts.begin();
    }

    template<class T>
    void Pool<T>::on_load_change() {
        std::lock_guard<std::mutex> guard(mutex);
        std::sort(
                ts.begin(),
                ts.end(),
                [](const auto &a,const auto &b) {
                    return a->current_load() < b->current_load();
                }
        );
    }

    template<class T>
    void Pool<T>::on_failure(const std::shared_ptr<T> t) {
        std::lock_guard<std::mutex> guard(mutex);
        std::remove_if(
                ts.begin(),
                ts.end(),
                [&](auto f) { return f == t; }
        );
    }
}