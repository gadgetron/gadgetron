#pragma once

#include <memory>
#include <istream>

#include "io/primitives.h"

namespace Gadgetron::Core {

    class StorageSpace {
    public:

        class StreamProvider {
        public:
            virtual ~StreamProvider() = default;

            virtual std::unique_ptr<std::istream> istream(const std::string& key) = 0;
            virtual std::unique_ptr<std::ostream> ostream(const std::string& key) = 0;
        };

        explicit StorageSpace(std::unique_ptr<StreamProvider> provider);

        template<class T>
        T fetch(const std::string& key) {
            auto istream = provider->istream(key);
            return IO::read<T>(istream);
        }

        template<class T>
        void store(const std::string& key, const T& t) {
            auto ostream = provider->ostream(key);
            IO::write(ostream, t);
        }

    private:
        std::unique_ptr<StreamProvider> provider;
    };

    struct Storage {
        StorageSpace session, scanner, debug;
    };
}


