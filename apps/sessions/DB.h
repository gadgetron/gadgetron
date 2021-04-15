#pragma once

#include <rocksdb/db.h>
#include <rocksdb/merge_operator.h>
#include <filesystem>
#include "Types.h"
#include <boost/date_time.hpp>
#include <nlohmann/json.hpp>
#include <boost/hana/adapt_struct.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/keys.hpp>
#include <boost/hana/at_key.hpp>

namespace boost::posix_time {

    using json = nlohmann::json;
    void to_json(json& j, const ptime& time){
        j = to_iso_extended_string(time);
    }

    void from_json(const json& j, ptime& time ){
        time = from_iso_extended_string(j.get<std::string>());
    }

}


namespace Gadgetron::Sessions::DB {
    template<class T>
    std::enable_if_t<boost::hana::Struct<T>::value>
    to_json(nlohmann::json &j, const T &x) {
        namespace hana = boost::hana;
        hana::for_each(hana::keys(x), [&](auto name) {
            j[hana::to<char const*>(name)] = hana::at_key(x, name);
        });
    }
    template<class T>
    std::enable_if_t<boost::hana::Struct<T>::value>
    from_json(const nlohmann::json& j, T& x) {
        namespace hana = boost::hana;
        hana::for_each(hana::keys(x),[&](auto name){
            //Nicest line in the entire codebase.
            hana::at_key(x,name) = j.at(hana::to<char const*>(name)).template get<std::decay_t<decltype(hana::at_key(x,name))>>();
        });
    }

    using json = nlohmann::json;

    class DBError : public std::runtime_error {
    public:
        DBError(rocksdb::Status &status) : std::runtime_error(status.ToString()) {}
    };

    std::shared_ptr<rocksdb::DB> create_database(const std::filesystem::path &path) {
        rocksdb::Options options;
        options.create_if_missing = true;
        options.create_missing_column_families = true;
        rocksdb::DB *database;
        rocksdb::Status status = rocksdb::DB::Open(options, path.string(), &database);
        if (!status.ok()) throw DBError(status);
        return std::shared_ptr<rocksdb::DB>(database);
    };


    class JsonStore {
    public:
        Core::optional<json> get(std::string_view key) {
            auto result = std::string{};
            auto status = database->Get(rocksdb::ReadOptions(), handle, key, &result);
            if (status.ok()) return json::from_msgpack(result);
            if (status.IsNotFound()) return Core::none;
            throw DBError(status);
        }

        json at(std::string_view key) {
            auto result = std::string{};
            rocksdb::Status status = database->Get(rocksdb::ReadOptions(), handle, key, &result);
            if (status.ok()) return json::from_msgpack(result);
            if (status.IsNotFound()) throw std::out_of_range("Key no found");
            throw DBError(status);
        }


        void set(std::string_view key, const json &value) {
            std::string buffer;
            json::to_msgpack(value, buffer);
            auto status = database->Put(rocksdb::WriteOptions(), handle, key, buffer);
            if (!status.ok()) throw DBError(status);
        }

        void merge(std::string_view key, const json &value) {
            std::string buffer;
            json::to_msgpack(value, buffer);
            auto status = database->Merge(rocksdb::WriteOptions(), handle, key, buffer);
            if (!status.ok()) throw DBError(status);
        }

        template<class Collection>
        void delete_keys(Collection &&keys_to_be_deleted) {
            rocksdb::WriteBatch batch;
            for (const auto &key : keys_to_be_deleted) {
                batch.Delete(handle, key);
            }
            auto status = database->Write(rocksdb::WriteOptions(), &batch);
            if (!status.ok()) throw DBError(status);
        }

        void delete_key(std::string_view key) {
            auto status = database->Delete(rocksdb::WriteOptions(), handle, key);
            if (!status.ok() && !status.IsNotFound()) throw DBError(status);
        }

        ~JsonStore() {
            if (database && handle)
                database->DestroyColumnFamilyHandle(handle);
        }

    private:
        rocksdb::ColumnFamilyHandle *handle;
        std::shared_ptr<rocksdb::DB> database;

    };


    class JsonListPrepender : public rocksdb::AssociativeMergeOperator {
    public:
        ~JsonListPrepender() override = default;

        bool Merge(const rocksdb::Slice &key, const rocksdb::Slice *existing_value, const rocksdb::Slice &value,
                   std::string *new_value, rocksdb::Logger *logger) const override {

            json json_value = json::from_msgpack(value.ToStringView());

            if (existing_value) {
                json array = json::from_msgpack(existing_value->ToStringView());
                array.insert(array.begin(), json_value);
                json::to_msgpack(array, *new_value);
                return true;
            }
            json::to_msgpack(json_value, *new_value);
            return true;
        }
    };

    struct BlobMeta {
        std::string blob_id;
        boost::posix_time::ptime deletion_time;
    };

    struct PendingWrite {
        std::string blob_id;
        boost::posix_time::ptime transaction_expiration;
        BlobMeta meta;
    };
    struct PendingRead {
        std::string blob_id;
        boost::posix_time::ptime transaction_expiration;
    };

}
    BOOST_HANA_ADAPT_STRUCT(Gadgetron::Sessions::DB::BlobMeta, blob_id, deletion_time);
    BOOST_HANA_ADAPT_STRUCT(Gadgetron::Sessions::DB::PendingWrite, blob_id, transaction_expiration, meta);
    BOOST_HANA_ADAPT_STRUCT(Gadgetron::Sessions::DB::PendingRead, blob_id, transaction_expiration);


namespace Gadgetron::Sessions::DB {

    template<class T, bool WithMerge = false>
    class KeyValueStore {
    public:
        Core::optional<T> operator[](std::string_view key){
            auto j = store.get(key);
            if (j) return j->get<T>();
            return Core::none;
        }

        T at(std::string_view key){
            return store.at(key).get<T>();
        }

        void set(std::string_view key, const T& value){
            store.set(key,json(value));
        }
        void delete_key(std::string_view key){
            store.delete_key(key);
        }
        template<class COLLECTION>
                void delete_keys(COLLECTION&& collection){
            store.template delete_keys(collection);
        }

        template< class Dummy = void>
        std::enable_if_t<WithMerge,Dummy> merge(std::string_view key, const T& new_value){
            store.merge(key,new_value);
                    }


    private:
        JsonStore store;
    };


}






