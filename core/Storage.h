#pragma once

#include <memory>
#include <istream>
#include <ostream>

#include "io/ismrmrd_types.h"
#include "io/adapt_struct.h"
#include "io/primitives.h"
#include <boost/date_time/posix_time/posix_time_config.hpp>
#include <boost/iostreams/device/array.hpp>

namespace Gadgetron::Core {

    class StorageSpace {
    public:
        template<class T>
        class StorageList {
        public:
            T operator[](size_t index) {
                auto data = storageSpace.provider->fetch(keys.at(index));
                return IO::read<T>(*istream_from_data(data));
            }


            size_t size() { return keys.size(); }

        private:
            friend StorageSpace;

            StorageList( StorageSpace& storageSpace, std::vector<std::string> keys) : keys(std::move(keys)),
                                                                                    storageSpace(storageSpace) {}

            std::vector<std::string> keys;
            StorageSpace &storageSpace;

        };

        template<class T>
        friend
        class StorageList;


        class StreamProvider {
        public:
            virtual ~StreamProvider() = default;

            [[nodiscard]] virtual std::vector<std::string> content(const std::string &key) const = 0;

            [[nodiscard]] virtual std::vector<char> fetch(const std::string &key) const = 0;

            virtual void store(const std::string &key, const std::vector<char> &data,
                               boost::posix_time::time_duration duration) const = 0;
        };

        StorageSpace(std::shared_ptr<StreamProvider> provider, boost::posix_time::time_duration default_duration);
        StorageSpace() = default;

        template<class T>
        StorageList<T> fetch(const std::string &key) {
            return StorageList<T>(*this, provider->content(key));
        }

        template<class T>
        void store(const std::string &key, const T &value) const {
            this->store(key, value, default_duration);
        }

        template<class T>
        void store(const std::string &key, const T &value, boost::posix_time::time_duration duration) const {
            std::vector<char> data{};
            auto os = ostream_view(data);
            IO::write(*os,value);
            os->flush();
            provider->store(key, data, duration);
        }

    private:
        std::shared_ptr<StreamProvider> provider;
        boost::posix_time::time_duration default_duration;


        static std::unique_ptr<std::istream> istream_from_data(const std::vector<char> &data);

        static std::unique_ptr<std::ostream> ostream_view(std::vector<char> &data);

    };


    struct Storage {
        StorageSpace session, scanner, debug;
    };
}


