#pragma once

#include <chrono>
#include <memory>
#include <istream>
#include <ostream>
#include <map>
#include <optional>

#include "io/adapt_struct.h"
#include "io/ismrmrd_types.h"
#include "io/primitives.h"
#include <boost/date_time/posix_time/posix_time_config.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace Gadgetron::Storage {
    
struct StorageItemTags {
public:
    class Builder;

    std::string subject;
    std::optional<std::string> device;
    std::optional<std::string> session;
    std::optional<std::string> name;
    std::multimap<std::string, std::string> custom_tags;
};

class StorageItemTags::Builder {
public:
    Builder(std::string const& subject) { tags.subject = subject; }

    Builder& with_device(std::string const& device) {
        tags.device = device;
        return *this;
    }

    Builder& with_session(std::string const& session) {
        tags.session = session;
        return *this;
    }

    Builder& with_name(std::string const& name) {
        tags.name = name;
        return *this;
    }

    Builder& with_custom_tag(std::string const& tag_name, std::string const& tag_value) {
        tags.custom_tags.insert({tag_name, tag_value});
        return *this;
    }

    StorageItemTags build() { return tags; }

private:
    StorageItemTags tags;
};

struct StorageItem {
public:
    StorageItemTags tags;
    std::string location;
    std::string contentType;
    std::chrono::time_point<std::chrono::system_clock> lastModified;
    std::optional<std::chrono::time_point<std::chrono::system_clock>> expires;
    std::string data;
};

struct StorageItemList {
public:
    std::vector<StorageItem> items;
    bool complete;
    std::string continuation;
};

class StorageClient {
public:
    class Builder;

    StorageClient(std::string base_url) : base_url(base_url) {}

    StorageItemList list_items(StorageItemTags const& tags, size_t limit = 20);

    StorageItemList get_next_page_of_items(StorageItemList const& page);

    std::shared_ptr<std::istream> get_latest_item(StorageItemTags const& tags);

    std::shared_ptr<std::istream> get_item_by_url(const std::string& url);

    StorageItem store_item(StorageItemTags const& tags, std::istream& data, std::optional<std::chrono::seconds> time_to_live = {});

    std::optional<std::string> health_check();

private:
    std::string base_url;
};

class StreamProvider {
  public:
    virtual ~StreamProvider() = default;

    [[nodiscard]] virtual std::vector<std::string> content(const std::string& subject,
                                                           const std::string& key) const = 0;

    [[nodiscard]] virtual std::vector<char> fetch(const std::string& uuid) const = 0;

    virtual void store(const std::string& subject, const std::string& key, const std::vector<char>& data,
                       boost::posix_time::time_duration duration) = 0;
};

std::unique_ptr<std::istream> istream_from_data(const std::vector<char>& data);
std::unique_ptr<std::ostream> ostream_view(std::vector<char>& data);

template <class T> class StorageList {
  public:
    T operator[](size_t index) {
        auto data = provider->fetch(keys.at(index));
        return Core::IO::read<T>(*istream_from_data(data));
    }

    StorageList& operator=(StorageList&&) noexcept = default;

    size_t size() { return keys.size(); }

    bool empty() { return keys.empty(); }

    StorageList(std::shared_ptr<StreamProvider> provider, std::vector<std::string> keys)
        : keys(std::move(keys)), provider(std::move(provider)) {}

    auto begin() { return boost::make_transform_iterator(keys.begin(), iterator_transform()); }

    auto end() { return boost::make_transform_iterator(keys.end(), iterator_transform()); }

  private:
    std::function<T(const std::string&)> iterator_transform() {
        return [this](const std::string& key) -> T {
            auto data = provider->fetch(key);
            return Core::IO::read<T>(*istream_from_data(data));
        };
    }

    std::vector<std::string> keys;
    std::shared_ptr<StreamProvider> provider;
};

class GenericStorageSpace {
  public:
    GenericStorageSpace(std::shared_ptr<StreamProvider> provider, const Core::optional<std::string>& subject,
                        boost::posix_time::time_duration default_duration);

    GenericStorageSpace() = default;

    template <class T> void store(const std::string& key, const T& value) { this->store(key, value, default_duration); }

    template <class T> void store(const std::string& key, const T& value, boost::posix_time::time_duration duration) {
        if (!subject)
            throw std::runtime_error("Storage space is unavailable due to missing information in the ISMRMRD header");
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
} // namespace Gadgetron::Storage

namespace Gadgetron {

class StorageSpace : public Storage::GenericStorageSpace {
  public:
    using GenericStorageSpace::GenericStorageSpace;

    template <class T> Storage::StorageList<T> fetch(const std::string& key) const {
        if (this->subject)
            return Storage::StorageList<T>(this->provider, this->provider->content(*subject, key));
        return Storage::StorageList<T>(this->provider, {});
    }
};

class MeasurementSpace : public Storage::GenericStorageSpace {
  public:
    using GenericStorageSpace::GenericStorageSpace;

    template <class T> Storage::StorageList<T> fetch(const std::string& measurementID, const std::string& key) const {
        return Storage::StorageList<T>(this->provider, this->provider->content(measurementID, key));
    }
};

struct StorageSpaces {
    StorageSpace session, scanner;
    MeasurementSpace measurement;
};

} // namespace Gadgetron
