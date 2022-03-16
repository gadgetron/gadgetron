#pragma once

#include <chrono>
#include <istream>
#include <map>
#include <memory>
#include <optional>
#include <ostream>

#include "IsmrmrdContextVariables.h"
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

    StorageClient(std::string base_url) { this->base_url = base_url.erase(base_url.find_last_not_of("/") + 1); }

    StorageItemList list_items(StorageItemTags const& tags, size_t limit = 20);

    StorageItemList get_next_page_of_items(StorageItemList const& page);

    std::shared_ptr<std::istream> get_latest_item(StorageItemTags const& tags);

    std::shared_ptr<std::istream> get_item_by_url(std::string const& url);

    StorageItem store_item(StorageItemTags const& tags, std::istream& data,
                           std::optional<std::chrono::seconds> time_to_live = {});

    std::optional<std::string> health_check();

  private:
    std::string base_url;
};

class GenericStorageSpace {
  public:
    GenericStorageSpace(std::shared_ptr<StorageClient> client, IsmrmrdContextVariables context_vars,
                        std::chrono::seconds default_duration)
        : client(client), context_vars(context_vars), default_duration(default_duration) {}

    template <typename Rep, typename Period>
    GenericStorageSpace(std::shared_ptr<StorageClient> client, StorageItemTags::Builder tag_builder,
                        std::chrono::duration<Rep, Period> default_duration)
        : GenericStorageSpace(client, tag_builder, std::chrono::duration_cast<std::chrono::seconds>(default_duration)) {
    }

    virtual ~GenericStorageSpace() = default;

    template <class T> void store(const std::string& key, const T& value) { this->store(key, value, default_duration); }

    template <class T, typename Rep, typename Period>
    void store(const std::string& key, const T& value, std::chrono::duration<Rep, Period> duration) {
        auto tags = get_tag_builder(true).with_name(key).build();
        std::stringstream stream;
        Core::IO::write(stream, value);
        client->store_item(tags, stream, std::chrono::duration_cast<std::chrono::seconds>(duration));
    }

  protected:
    virtual StorageItemTags::Builder get_tag_builder(bool for_write) = 0;

    template <class T> std::optional<T> get_latest(StorageItemTags const& tags) const {
        auto data = client->get_latest_item(tags);
        if (data) {
            return Core::IO::read<T>(*data);
        }

        return {};
    }

    std::shared_ptr<StorageClient> client;
    IsmrmrdContextVariables context_vars;
    std::chrono::seconds default_duration;
};
} // namespace Gadgetron::Storage

namespace Gadgetron {

using namespace Gadgetron::Storage;

class IncompleteStorageContextException : std::exception {
  public:
    IncompleteStorageContextException(std::string what) : what_(what) {}

    const char* what() const noexcept override { return what_.c_str(); }

  protected:
    std::string what_;
};

class StorageSpaceWithDefaultRead : public GenericStorageSpace {
  public:
    using GenericStorageSpace::GenericStorageSpace;

    template <class T> std::optional<T> get_latest(const std::string& key) {
        try {
            auto tags = get_tag_builder(false).with_name(key).build();
            return GenericStorageSpace::get_latest<T>(tags);
        } catch (IncompleteStorageContextException const&) {
            return {};
        }
    }
};

class SessionSpace : public StorageSpaceWithDefaultRead {
  public:
    using StorageSpaceWithDefaultRead::StorageSpaceWithDefaultRead;

  protected:
    StorageItemTags::Builder get_tag_builder(bool for_write) override {
        if (context_vars.subject_id().empty()) {
            throw IncompleteStorageContextException(
                "Storage space is unavailable due to missing information in the ISMRMRD header");
        }

        auto builder = StorageItemTags::Builder(context_vars.subject_id());
        if (!context_vars.device_id().empty()) {
            builder.with_device(context_vars.device_id());
        }

        if (!context_vars.session_id().empty()) {
            builder.with_session(context_vars.session_id());
        }

        return builder;
    }
};

class ScannerSpace : public StorageSpaceWithDefaultRead {
  public:
    using StorageSpaceWithDefaultRead::StorageSpaceWithDefaultRead;

  protected:
    StorageItemTags::Builder get_tag_builder(bool for_write) override {
        if (context_vars.device_id().empty()) {
            throw IncompleteStorageContextException(
                "Storage space is unavailable due to missing information in the ISMRMRD header");
        }

        return StorageItemTags::Builder("$null").with_device(context_vars.device_id());
    }
};

class MeasurementSpace : public GenericStorageSpace {
  public:
    using GenericStorageSpace::GenericStorageSpace;

    template <class T> std::optional<T> get_latest(const std::string& measurement_id, const std::string& key) {
        try {
            auto tags = get_tag_builder(false).with_name(key).with_custom_tag("measurement", measurement_id).build();
            return GenericStorageSpace::get_latest<T>(tags);
        } catch (IncompleteStorageContextException const&) {
            return {};
        }
    }

  protected:
    StorageItemTags::Builder get_tag_builder(bool for_write) override {
        auto subject_id = context_vars.subject_id().empty() ? "$null" : context_vars.subject_id();
        auto builder = StorageItemTags::Builder(subject_id);

        if (context_vars.measurement_id().empty()) {
            if (for_write) {
                throw IncompleteStorageContextException(
                    "Storage space is unavailable due to missing information in the ISMRMRD header");
            }
            // reads will use the dependency measurement id which is not taken from the context
        } else {
            builder.with_custom_tag("measurement", context_vars.measurement_id());
        }

        if (!context_vars.device_id().empty()) {
            builder.with_device(context_vars.device_id());
        }

        if (!context_vars.session_id().empty()) {
            builder.with_session(context_vars.session_id());
        }

        return builder;
    }
};

struct StorageSpaces {
    std::shared_ptr<SessionSpace> session;
    std::shared_ptr<ScannerSpace> scanner;
    std::shared_ptr<MeasurementSpace> measurement;
};

} // namespace Gadgetron
