#pragma once

#include <memory>
#include <istream>

#include "io/ismrmrd_types.h"
#include "io/adapt_struct.h"
#include "io/primitives.h"

namespace Gadgetron::Core {



    class StorageSpace {
    public:

        class StreamProvider {
        public:
            virtual ~StreamProvider() = default;

            virtual std::unique_ptr<std::istream> fetch(const std::string& key) const = 0;
            virtual void store(const std::string& key, std::istream& datastream) const = 0;
            virtual bool contains(const std::string& key) const = 0;
        };

        explicit StorageSpace(std::shared_ptr<StreamProvider> provider);

        StorageSpace() {}

        template<class T>
        T fetch(const std::string& key) const {
            auto istream = provider->fetch(key);
            return IO::read<T>(*istream);
        }

        template<class T>
        void store(const std::string& key, const T& t) const {
            std::stringstream stream;
            IO::write(stream, t);
            provider->store(key,stream);

        }

        bool contains(const std::string& key) const {
            return provider->contains(key);
        }

    private:
        std::shared_ptr<StreamProvider> provider;
    };

    struct Storage {
        StorageSpace session, noise, debug;
    };
}


