#pragma once

#include "JSONStore.h"
#include "DBError.h"
#include <rocksdb/db.h>
#include <rocksdb/merge_operator.h>
#include <boost/filesystem.hpp>
#include <boost/date_time.hpp>
#include <nlohmann/json.hpp>
#include <boost/hana/adapt_struct.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/keys.hpp>
#include <boost/hana/at_key.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/string_generator.hpp>

namespace boost::posix_time {

    using json = nlohmann::json;

    inline static void to_json(json &j, const ptime &time) {
        j = to_iso_extended_string(time);
    }

    inline static void from_json(const json &j, ptime &time) {
        time = from_iso_extended_string(j.get<std::string>());
    }

    inline static void to_json(json &j, const time_duration &duration) {
        j = to_simple_string(duration);
    }

    inline void from_json(const json &j, time_duration &duration) {
        duration = duration_from_string(j.get<std::string>());
    }

}
namespace boost::uuids {

    using json = nlohmann::json;

    inline static void to_json(json &j, const uuid &id) {
        j = to_string(id);
    }

    inline static void from_json(const json &j, uuid &id) {
        boost::uuids::string_generator gen;
        id = gen(j.get<std::string>());
    }
}

namespace nlohmann {
    template<class T>
    std::enable_if_t<boost::hana::Struct<T>::value>
    to_json(nlohmann::json &j, const T &x) {
        namespace hana = boost::hana;
        hana::for_each(hana::keys(x), [&](auto name) {
            j[hana::to<char const *>(name)] = hana::at_key(x, name);
        });
    }

    template<class T>
    std::enable_if_t<boost::hana::Struct<T>::value>
    from_json(const nlohmann::json &j, T &x) {
        namespace hana = boost::hana;
        hana::for_each(hana::keys(x), [&](auto name) {
            //Nicest line in the entire codebase.
            hana::at_key(x, name) = j.at(hana::to<char const *>(name)).template get<std::decay_t<decltype(hana::at_key(
                    x, name))>>();
        });
    }
}


namespace Gadgetron::Storage::DB {


    using json = nlohmann::json;



    struct DataBaseFamilies {
        std::shared_ptr<rocksdb::DB> database;
        std::map<std::string, rocksdb::ColumnFamilyHandle *> families;
    };

    inline DataBaseFamilies create_database(const boost::filesystem::path &path,
                                            const std::vector<rocksdb::ColumnFamilyDescriptor> &column_families) {
        rocksdb::DBOptions options;
        options.create_if_missing = true;
        options.create_missing_column_families = true;

        options.db_log_dir = "/tmp/dblog";
        options.OptimizeForSmallDb();
        rocksdb::DB *database;
        std::vector<rocksdb::ColumnFamilyHandle *> handles;


        rocksdb::Status status = rocksdb::DB::Open(options, path.string(), column_families, &handles, &database);
        if (!status.ok()) throw DBError(status);

        return DataBaseFamilies{std::shared_ptr<rocksdb::DB>(database), ranges::zip_view(
                column_families | ranges::views::transform([](const auto &cf) { return cf.name; }), handles) |
                                                                        ranges::to<std::map>()};
    };


    struct BlobMeta {
        boost::uuids::uuid blob_id;
        boost::posix_time::ptime creation_time;
        boost::posix_time::ptime deletion_time;

        bool operator==(const BlobMeta &other) const {
            auto compare_timestamps = [](auto &time1, auto &time2) {
                return boost::posix_time::to_iso_extended_string(time1) ==
                       boost::posix_time::to_iso_extended_string(time2);
            };

            return (blob_id == other.blob_id) && compare_timestamps(creation_time, other.creation_time) &&
                   compare_timestamps(deletion_time, other.deletion_time);
        }
    };

    struct PendingWrite {
        std::string key;
        boost::posix_time::ptime transaction_expiration;
        BlobMeta meta;
    };


}
BOOST_HANA_ADAPT_STRUCT(Gadgetron::Storage::DB::BlobMeta, blob_id, creation_time, deletion_time);
BOOST_HANA_ADAPT_STRUCT(Gadgetron::Storage::DB::PendingWrite, key, transaction_expiration, meta);


namespace Gadgetron::Storage::DB {

    struct Database {
        static std::shared_ptr<Database> make_db(const boost::filesystem::path &database_path) {
            auto cf_descriptors = std::vector<rocksdb::ColumnFamilyDescriptor>{{"Info", rocksdb::ColumnFamilyOptions()},
                                                                               {"PendingWrites", rocksdb::ColumnFamilyOptions()},
                                                                               {"Blobs", rocksdb::ColumnFamilyOptions()},
                                                                               {"default", rocksdb::ColumnFamilyOptions()}};
            auto families = create_database(database_path, cf_descriptors);
            return std::make_shared<Database>(families);
        }
        Database(DataBaseFamilies &families) : db_info(families.database, families.families.at("Info")),
                                               pending_writes(families.database, families.families.at("PendingWrites")),
                                               blobs(families.database, families.families.at("Blobs")) {

        }

        JSONStore db_info;
        ValueStore<PendingWrite> pending_writes;
        ListStore<BlobMeta> blobs;

    };

}






