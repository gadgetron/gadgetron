#pragma once

#include <rocksdb/db.h>
#include <filesystem>
#include "Types.h"
#include <boost/date_time.hpp>
#include <nlohmann/json.hpp>
namespace Gadgetron::Sessions::DB {
    namespace json = nlohmann::json;
    class DBError : public std::runtime_error {
    public:
        DBError(rocksdb::Status &status) : std::runtime_error(status.ToString()) {}
    };

    std::shared_ptr<rocksdb::DB> create_database(const std::filesystem::path &path) {
        rocksdb::Options options;
        options.create_if_missing = true;
        options.create_missing_column_families = true;
        std::shared_ptr<rocksdb::DB> database;
        rocksdb::Status status = rocksdb::DB::Open(options, path.string(), &database.get());
        if (!status.ok()) throw DBError(status);
        return database;
    };


    class ColumnFamily {
    public:
        Core::Optional<std::string> get(const rocksdb::Slice& key){
            auto result = std::string{};
            auto status = database->Get(rocksdb::ReadOptions(),handle,key, &result);
            if (status.ok()) return  result;
            if (status.IsNotFound()) return Core::none;
            throw DBError(status);
        }
        std::string at(const rocksdb::Slice& key){
            auto result = std::string{};
            rocksdb::Status status = database->Get(rocksdb::ReadOptions(),handle,key, &result);
            if (status.ok()) return  result;
            if (status.IsNotFound()) throw std::out_of_range("Key no found");
            throw DBError(status);
        }


        void set(const rocksdb::Slice& key , const rocksdb::Slice& value){
            auto status = database->Put(rocksdb::WriteOptions(),handle,key,value);
            if (!status.ok()) throw DBError(status);
        }

        void merge(const rocksdb::Slice& key, const rocksdb::Slice& value ){
            auto status = database->Merge(rocksdb::WriteOptions(),handle,key,value);
            if (!status.ok()) throw DBError(status);
        }

        ~ColumnFamily(){
            if (database && handle)
                database->DestroyColumnFamilyHandle(handle);
        }

    private:
        rocksdb::ColumnFamilyHandle* handle;
        std::shared_ptr<rocksdb::DB> database;

    };

    struct BlobMeta {
        std::string blob_id;
        boost::posix_time::ptime expiration_time;
    };

    class BlobDB {
    public:
        Core::optional<BlobMeta> get(std::string_view key ){
            using namespace rocksdb
            auto result = column->get(Slice(key));
            if (!result) return Core::none;
           auto j = json::parse(result.get());
           return j.get<BlobMeta>();


        }

        void set(std::string_view key, const BlobMeta& data){
            using namespace rocksdb;

            json::json j = data;
            std::string serialized = j.dump();
            column->set(Slice(key),Slice(serialized));
        }
    private:
        std::unique_ptr<ColumnFamily> column;
    };

    struct PendingTransaction{
        std::string key;
        std::string blob_id;
    };



}