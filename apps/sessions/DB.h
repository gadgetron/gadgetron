#pragma once

#include "JSONStore.h"
#include "DBError.h"
#include <rocksdb/db.h>
#include <rocksdb/merge_operator.h>
#include <filesystem>
#include <boost/date_time.hpp>
#include <nlohmann/json.hpp>
#include <boost/hana/adapt_struct.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/keys.hpp>
#include <boost/hana/at_key.hpp>

namespace boost::posix_time {

    using json = nlohmann::json;

    void to_json(json &j, const ptime &time) {
        j = to_iso_extended_string(time);
    }

    void from_json(const json &j, ptime &time) {
        time = from_iso_extended_string(j.get<std::string>());
    }

}


namespace Gadgetron::Sessions::DB {
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

    using json = nlohmann::json;



    struct DataBaseFamilies {
        std::shared_ptr<rocksdb::DB> database;
        std::map<std::string,rocksdb::ColumnFamilyHandle*> families;
    };

    DataBaseFamilies create_database(const std::filesystem::path &path, const std::map<std::string,rocksdb::ColumnFamilyDescriptor>& column_families) {
        rocksdb::DBOptions options;
        options.create_if_missing = true;
        options.create_missing_column_families = true;
        rocksdb::DB *database;
        auto column_descriptors = column_families | ranges::views::values | ranges::to<std::vector>();
        std::vector<rocksdb::ColumnFamilyHandle*> handles;


        rocksdb::Status status = rocksdb::DB::Open(options, path.string(), column_descriptors,&handles, &database);
        if (!status.ok()) throw DBError(status);


        return DataBaseFamilies{std::shared_ptr<rocksdb::DB>(database), ranges::zip_view(column_families | ranges::views::keys , handles) | ranges::to<std::map>()};
    };





    struct BlobMeta {
        std::string blob_id;
        boost::posix_time::ptime creation_time;
        boost::posix_time::ptime deletion_time;
    };

    struct PendingWrite {
        std::string blob_id;
        boost::posix_time::ptime transaction_expiration;
        BlobMeta meta;
    };


}
BOOST_HANA_ADAPT_STRUCT(Gadgetron::Sessions::DB::BlobMeta, blob_id,creation_time, deletion_time);
BOOST_HANA_ADAPT_STRUCT(Gadgetron::Sessions::DB::PendingWrite, blob_id, transaction_expiration, meta);


namespace Gadgetron::Sessions::DB {



}






