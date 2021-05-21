#pragma once


#include "Types.h"
#include <rocksdb/db.h>
#include <rocksdb/iterator.h>
#include <nlohmann/json.hpp>
#include <boost/iterator.hpp>
#include "DBError.h"
#include <range/v3/view.hpp>
#include <range/v3/range_concepts.hpp>
#include <range/v3/algorithm/any_of.hpp>


namespace Gadgetron::Storage::DB {
    using json = nlohmann::json;


    template<class T>
    class ValueStore;

    template<class T>
    class Iterator {

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = std::pair<std::string, T>;
        using difference_type = int; //We really can't take the difference between these guys
        using pointer = const value_type *;
        using reference = value_type;

        Iterator() = default;

        Iterator(const Iterator &) = default;

        bool operator==(const Iterator &other) const {
            if (other.it == this->it) return true;
            if (other.it) return false;
            //Empty iterator is used as sentinel
            if (it->Valid()) {
                return false;
            }
            if (!it->status().ok()) {
                throw DBError(it->status());
            }
            return true;
        }

        bool operator!=(const Iterator &e) const {
            return !(*this == e);
        }

        reference operator*() const {
            return {it->key().ToString(), json::from_msgpack(it->value().ToStringView()).get<T>()};
        }

        void operator++(int) {
            it->Next();
        }

        Iterator &operator++() {
            it->Next();
            return *this;
        }

        Iterator &operator--() {
            it->Prev();
            return *this;
        }

    private:
        explicit Iterator(rocksdb::Iterator *itr_ptr) : it(itr_ptr) { it->SeekToFirst(); };

        Iterator(rocksdb::Iterator *itr_ptr, std::string_view from) : it(itr_ptr) { it->Seek(from); };
        std::shared_ptr<rocksdb::Iterator> it;

        friend class JSONStore;

        friend class ValueStore<T>;
    };


    class JSONStore {
    public:

        JSONStore(const JSONStore &) = delete;

        JSONStore(JSONStore &&) = default;

        JSONStore(std::shared_ptr<rocksdb::DB> database, rocksdb::ColumnFamilyHandle *handle) : database(
                std::move(database)), handle(handle) {}

        JSONStore &operator=(JSONStore &&other) noexcept {
            handle = other.handle;
            database = std::move(other.database);
            other.handle = nullptr;
            return *this;

        }

        Core::optional<json> operator[](std::string_view key) const {
            auto result = std::string{};
            auto status = database->Get(rocksdb::ReadOptions(), handle, key, &result);
            if (status.ok()) return json::from_msgpack(result);
            if (status.IsNotFound()) return Core::none;
            throw DBError(status);
        }

        json at(std::string_view key) const {
            auto result = std::string{};
            rocksdb::Status status = database->Get(rocksdb::ReadOptions(), handle, key, &result);
            if (status.ok()) return json::from_msgpack(result);
            if (status.IsNotFound()) throw std::out_of_range("Key no found");
            throw DBError(status);
        }

        void set(std::string_view key, const json &value) {
            std::string str_val;
            json::to_msgpack(value, str_val);
            auto options = rocksdb::WriteOptions();
            auto status = database->Put(options, handle, key, str_val);
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

        template<class RANGE>
        void update(RANGE &&key_values) {
            rocksdb::WriteBatch batch;
            for (const auto&[key, value] : key_values) {
                std::string buffer;
                json::to_msgpack(value, buffer);
                batch.Put(handle, key, buffer);
            }
            database->Write(rocksdb::WriteOptions(), &batch);
        }

        ~JSONStore() {
            if (database && handle)
                database->DestroyColumnFamilyHandle(handle);
        }


        Iterator<json> begin() {
            auto it = database->NewIterator(rocksdb::ReadOptions(), handle);
            return Iterator<json>(it);
        }

        Iterator<json> begin(std::string_view from) {
            auto it = database->NewIterator(rocksdb::ReadOptions(), handle);
            return Iterator<json>(it, from);
        }

        Iterator<json> end() {
            return {};
        }


        auto create_db_iterator() {
            return database->NewIterator(rocksdb::ReadOptions(), handle);
        }


    private:

        mutable rocksdb::ColumnFamilyHandle *handle = nullptr;
        mutable std::shared_ptr<rocksdb::DB> database = nullptr;

    };


    template<class T>
    class ValueStore {

    public:
        ValueStore(std::shared_ptr<rocksdb::DB> ptr, rocksdb::ColumnFamilyHandle *handle) : store(std::move(ptr),
                                                                                                  handle) {}

        ValueStore() = default;

        ValueStore(const ValueStore &) = delete;

        ValueStore &operator=(ValueStore &&other) noexcept {
            store = std::move(other.store);
            return *this;
        }

        Core::optional<T> operator[](std::string_view key) {
            if (auto j = store[key]) {
                return j->get<T>();
            }
            return {};
        }

        auto begin() {
            return Iterator<T>(store.create_db_iterator());
        }

        auto end() {
            return Iterator<T>();
        }

        void set(std::string_view key, const T &val) {
            store.set(key, json(val));
        }

        template<class COLLECTION>
        void delete_keys(COLLECTION &&collection) {
            store.template delete_keys(std::forward<COLLECTION>(collection));
        }

        void delete_key(std::string_view str) {
            store.delete_key(str);
        }

        template<class RANGE>
        void update(RANGE &&key_values) {

            store.template update(key_values | ranges::views::transform([](const auto &key_value) {
                const auto&[key, value] = key_value;
                return std::make_pair(key, json(value));
            }));
        }


    private:
        JSONStore store;
    };


    template<class T>
    class ListStore {
    public:
        ListStore() = default;

        ListStore(const ListStore &) = delete;

        ListStore(std::shared_ptr<rocksdb::DB> ptr, rocksdb::ColumnFamilyHandle *handle) : store(std::move(ptr),
                                                                                                 handle) {}

        ListStore &operator=(ListStore &&other) noexcept {
            store = std::move(other.store);
            return *this;
        }

        void push_back(std::string_view key, const T &val) {
            auto list = store[key];
            if (list) {
                list->push_back(val);
                store.set(key, *list);
                return;
            }
            store.set(key, std::vector<T>{val});
        }

        auto begin() {
            return store.begin();
        }

        auto end() {
            return store.end();
        }


        template<class RANGE>
        void update(RANGE &&key_value_range) {
            store.template update(std::forward<RANGE>(key_value_range));
        }


        std::vector<T> operator[](std::string_view key) {
            auto list = store[key];
            if (list) return std::move(*list);
            return {};
        }

    private:
        ValueStore<std::vector<T>> store;


    };
}
