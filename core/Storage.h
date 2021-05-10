#pragma once

#include <memory>
#include <istream>
#include <ostream>

#include "io/ismrmrd_types.h"
#include "io/adapt_struct.h"
#include "io/primitives.h"
#include <boost/date_time/posix_time/posix_time_config.hpp>
#include <boost/iostreams/device/array.hpp>

namespace Gadgetron::Storage {
    class StreamProvider {
    public:
        virtual ~StreamProvider() = default;

        [[nodiscard]] virtual std::vector<std::string>
        content(const std::string &subject, const std::string &key) const = 0;

        [[nodiscard]] virtual std::vector<char> fetch(const std::string &uuid) const = 0;

        virtual void store(const std::string &subject, const std::string &key, const std::vector<char> &data,
                           boost::posix_time::time_duration duration) = 0;
    };

    std::unique_ptr<std::istream> istream_from_data(const std::vector<char> &data);

    std::unique_ptr<std::ostream> ostream_view(std::vector<char> &data);

    template<class T>
    class StorageList {
    public:
        T operator[](size_t index) {
            auto data = provider->fetch(keys.at(index));
            return Core::IO::read<T>(*istream_from_data(data));
        }

        StorageList &operator=(StorageList &&) noexcept = default;

        size_t size() { return keys.size(); }

        bool empty(){return keys.empty();}

        StorageList(std::shared_ptr<StreamProvider> provider, std::vector<std::string> keys) : keys(std::move(keys)),
                                                                                               provider(std::move(
                                                                                                       provider)) {}

    private:
        std::vector<std::string> keys;
        std::shared_ptr<StreamProvider> provider;

    };

    class GenericStorageSpace {
    public:

        GenericStorageSpace(std::shared_ptr<StreamProvider> provider, const Core::optional<std::string> &subject,
                            boost::posix_time::time_duration default_duration);

        GenericStorageSpace() = default;

        template<class T>
        void store(const std::string &key, const T &value) {
            this->store(key, value, default_duration);
        }

        template<class T>
        void store(const std::string &key, const T &value, boost::posix_time::time_duration duration) {
            if (!subject)
                throw std::runtime_error(
                        "Storage space is unavailable due to missing information in the ISMRMRD header");
            std::vector<char> data{};
            auto os = ostream_view(data);
            Core::IO::write(*os, value);
            os->flush();
            provider->store(*subject, key, data, duration);
        }

    protected:
        Core::optional<std::string> subject;
        std::shared_ptr<StreamProvider> provider;
        boost::posix_time::time_duration default_duration;
    };
}

namespace Gadgetron {

class StorageSpace : public Storage::GenericStorageSpace {
    public:
        using GenericStorageSpace::GenericStorageSpace;

        template<class T>
        Storage::StorageList<T> fetch(const std::string &key) const {
            if (this->subject)
                return Storage::StorageList<T>(this->provider,this->provider->content(*subject, key));
            return Storage::StorageList<T>(this->provider,{});
        }
    };

    class MeasurementSpace : public Storage::GenericStorageSpace {
    public:
        using GenericStorageSpace::GenericStorageSpace;

        template<class T>
        Storage::StorageList<T> fetch(const std::string &measurementID, const std::string &key) const {
            return Storage::StorageList<T>(this->provider,this->provider->content(measurementID, key));
        }
    };


    struct StorageSpaces {
        StorageSpace session, scanner, debug;
        MeasurementSpace measurment;
    };

}


