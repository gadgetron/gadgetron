#pragma once

#include <vector>
#include <memory>

namespace Gadgetron::Server::Connection::Stream {

    class CloseGuard {
    public:

        template<class ...ARGS>
        explicit CloseGuard(ARGS&... args) : closers{std::make_shared<Closeable<ARGS>>(args)...} {}

        template<class T>
        explicit CloseGuard(std::vector<T> &args) {
            std::transform(
                    args.begin(), args.end(),
                    std::back_inserter(closers),
                    [](auto t) { return std::make_shared<Closeable>(t); }
            );
        }

        template<class K, class V>
        explicit CloseGuard(std::map<K, V> &args) {
            std::transform(
                    args.begin(), args.end(),
                    std::back_inserter(closers),
                    [](auto pair) { return std::make_shared<Closeable>(pair.second); }
            );
        }

        ~CloseGuard() {
            for (auto &closer : closers) closer->close();
        }

    private:

        class ICloseable {
        public:
            virtual void close() = 0;
            virtual ~ICloseable() = default;
        };

        template<class T>
        class Closeable : public ICloseable {
        public:
            explicit Closeable(T &t) : t(t) {}
            void close() override { t.close(); }
        private:
            T &t;
        };

        std::vector<std::shared_ptr<ICloseable>> closers;
    };
}