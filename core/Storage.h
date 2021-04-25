#pragma once

#include <memory>
#include <istream>
#include <ostream>

#include "io/ismrmrd_types.h"
#include "io/adapt_struct.h"
#include "io/primitives.h"

namespace Gadgetron::Core {




    class StorageSpace {
    public:
        template<class T>
        class StorageList {

            T operator[](size_t index) {
                return storageSpace.fetch<T>(keys.at(index));
            }


            size_t size(){ return keys.size();}
        private:
            friend StorageSpace;
            StorageList(std::vector<std::string> keys, StorageSpace storageSpace): keys(std::move(keys)), storageSpace(storageSpace) {}
            std::vector<std::string> keys;
            StorageSpace& storageSpace;

        };

        template<class T>
        friend class StorageList;


        class StreamProvider {
        public:
            virtual ~StreamProvider() = default;
            virtual std::vector<std::string> content(const std::string& key) const = 0;
            virtual std::unique_ptr<std::istream> fetch(const std::string& key) const = 0;
            virtual std::unique_ptr<std::ostream> store(const std::string& key) const = 0;
        };

        explicit StorageSpace(std::shared_ptr<StreamProvider> provider);

        StorageSpace() {}

        template<class T>
        StorageList<T> fetch(const std::string& key){
            return StorageList<T>(*this,provider->content(key));
        }


        template<class T>
        void store(const std::string& key, const T& t) const {
            auto stream = provider->store(key);
            IO::write(*stream, t);
        }


    protected:

        template<class T>
        T fetch_content(const std::string& key) const {
            auto istream = provider->fetch(key);
            return IO::read<T>(*istream);
        }

    private:
        std::shared_ptr<StreamProvider> provider;
    };


    struct Storage {
        StorageSpace session, noise, debug;
    };
}


