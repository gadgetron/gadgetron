#pragma once


#include "Types.h"
#include <rocksdb/db.h>
#include <rocksdb/iterator.h>
#include <nlohmann/json.hpp>
#include <boost/iterator.hpp>
#include "DBError.h"
#include <range/v3/view.hpp>
#include <range/v3/range_concepts.hpp>



namespace Gadgetron::Storage::DB {
    using json = nlohmann::json;



    class JSONStore {
    public:
        JSONStore() = default;
        JSONStore(const JSONStore&) = delete;
        JSONStore(JSONStore&& ) = default;
        JSONStore(std::shared_ptr<rocksdb::DB> database, rocksdb::ColumnFamilyHandle* handle) : database(std::move(database)), handle(handle ){}

        JSONStore& operator=(JSONStore&& other) noexcept{
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

        ~JSONStore() {
            if (database && handle)
                database->DestroyColumnFamilyHandle(handle);
        }

        class Iterator;

        class ReverseIterator;

        class EndPoint {
        public:
            EndPoint(const EndPoint&) = default;
            EndPoint() : sentinel(Core::none) {}

            EndPoint(std::string_view end_string) : sentinel(std::in_place, end_string) {}

        private:
            Core::optional<std::string> sentinel;
            friend Iterator;
            friend ReverseIterator;
        };

        class Iterator {

        public:
            using iterator_category = std::input_iterator_tag;
            using value_type = std::pair<std::string, json>;
            using difference_type = int; //We really can't take the difference between these guys
            using pointer = const value_type *;
            using reference = value_type;

            Iterator() = default;
            Iterator(const Iterator &) = default;

            bool operator==(const EndPoint &endpoint) const {
                if (it->Valid()) {
                    if (endpoint.sentinel) {
                        return it->key().ToString() > endpoint.sentinel.value();
                    }
                    return false;
                }
                if (!it->status().ok()){
                    throw DBError(it->status());
                }

                return true;
            }

            bool operator!=(const EndPoint &e) const {
                return !(*this == e);
            }

            reference operator*() {
                return {it->key().ToString(), json::from_msgpack(it->value().ToStringView())};
            }

            Iterator &operator++() {
                it->Next();
                return *this;
            }

            Iterator operator++(int) { it->Next(); return *this;};

            Iterator &operator--() {
                it->Prev();
                return *this;
            }

            void operator--(int) { it->Prev();}

        private:
            explicit Iterator(rocksdb::Iterator *itr_ptr) : it(itr_ptr) { it->SeekToFirst(); };

            Iterator(rocksdb::Iterator *itr_ptr, std::string_view from) : it(itr_ptr) { it->Seek(from); };
            std::shared_ptr<rocksdb::Iterator> it;
            friend JSONStore;
        };


        class ReverseIterator {

        public:
            using iterator_category = std::input_iterator_tag;
            using value_type = std::pair<std::string, json>;
            using difference_type = void; //We really can't take the difference between these guys
            using pointer = const value_type *;
            using reference = value_type;

            ReverseIterator() = default;
            ReverseIterator(const ReverseIterator &) = default;

            bool operator==(const EndPoint &endpoint) const {
                if (it->Valid()) {
                    if (endpoint.sentinel) {
                        return it->key().ToString() < endpoint.sentinel.value();
                    }
                    return false;
                }

                return true;
            }

            bool operator!=(const EndPoint &e) const {
                return !(*this == e);
            }

            reference operator*() {
                return {it->key().ToString(), json::from_msgpack(it->value().ToStringView())};
            }

            ReverseIterator &operator++() {
                it->Prev();
                return *this;
            }

            void operator++(int) {
                it->Prev();
            }

            void operator--(int){
                it->Next();
            }

            ReverseIterator &operator--() {
                it->Next();
                return *this;
            }

        private:
            explicit ReverseIterator(rocksdb::Iterator *itr_ptr) : it(itr_ptr) { it->SeekToLast(); };

            ReverseIterator(rocksdb::Iterator *itr_ptr, std::string_view from) : it(itr_ptr) { it->SeekForPrev(from); };
            std::shared_ptr<rocksdb::Iterator> it;
            friend JSONStore;
        };


        Iterator begin() {
            auto it = database->NewIterator(rocksdb::ReadOptions(), handle);
            return Iterator(it);
        }

        Iterator begin(std::string_view from) {
            auto it = database->NewIterator(rocksdb::ReadOptions(), handle);
            return Iterator(it, from);
        }

        EndPoint end() {
            return EndPoint{};
        }

        EndPoint end(std::string_view to) {
            return EndPoint(to);
        }


        ReverseIterator rbegin() {
            auto it = database->NewIterator(rocksdb::ReadOptions(), handle);
            return ReverseIterator(it);
        }

        ReverseIterator rbegin(std::string_view from) {
            auto it = database->NewIterator(rocksdb::ReadOptions(), handle);
            return ReverseIterator(it, from);
        }

        EndPoint rend() {
            return EndPoint{};
        }

        EndPoint rend(std::string_view to) {
            return EndPoint(to);
        }
        friend class Range;

        struct Range {
            std::string lower_key, upper_key; JSONStore& store;
            Iterator begin() {
                auto options = rocksdb::ReadOptions();
                auto it = store.database->NewIterator(options,store.handle);
                return Iterator(it, lower_key);
            }

            EndPoint end() { return EndPoint(upper_key); };
        };

        friend class RevRange;
        struct RevRange {

            std::string lower_key, upper_key; JSONStore& store;
            ReverseIterator begin() {
                auto options = rocksdb::ReadOptions();
                auto it = store.database->NewIterator(options,store.handle);
                return ReverseIterator(it,upper_key);
            }
            EndPoint end() { return EndPoint(lower_key); };
        };

       /**
         *
         * @param from Inclusive
         * @param to Inclusive
         * @return
         */
        auto range(std::string_view from, std::string_view to) {
            assert(from<to);
            return Range{std::string(from), std::string(to),*this};
        }

        /**
        *
        * @param from Inclusive
        * @param to Inclusive
        * @return
        */
        auto reverse_range(std::string_view from, std::string_view to) {
            assert(from < to);
            return RevRange{std::string(from),std::string(to),*this};

        }

    public:

        mutable rocksdb::ColumnFamilyHandle *handle = nullptr;
        mutable std::shared_ptr<rocksdb::DB> database = nullptr;

    };


    template<class T>
    class ValueStore {

    public:
        ValueStore(std::shared_ptr<rocksdb::DB> ptr,rocksdb::ColumnFamilyHandle* handle) : store(std::move(ptr),handle) {}
        ValueStore() = default;
        ValueStore(const ValueStore&) = delete;

        ValueStore& operator=(ValueStore&& other) noexcept {
            store = std::move(other.store);
            return *this;
        }

        Core::optional<T> operator[](std::string_view key){
            if (auto j = store[key]){
                return j->get<T>();
            }
            return {};
        }

        void set(std::string_view key, const T& val){
            store.set(key,json(val));
        }

        template<class COLLECTION>
        void delete_keys(COLLECTION&& collection){
            store.template delete_keys(std::forward<COLLECTION>(collection));
        }

        void delete_key(std::string_view str){
            store.delete_key(str);
        }

    private:
        JSONStore store;
    };


    template<class T>
    class ListStore {
    public:
        ListStore() = default;
        ListStore(const ListStore&) = delete;
        ListStore(std::shared_ptr<rocksdb::DB> ptr,rocksdb::ColumnFamilyHandle* handle) : store(std::move(ptr),handle) {}

        ListStore& operator=(ListStore&& other) noexcept {
            store = std::move(other.store);
            return *this;
        }
        void push_back(std::string_view key, const T &val) {
            auto list = store[key];
            if (list){
                list->push_back(val);
                store.set(key,*list);
                return;
            }
            store.set(key,std::vector<T>{val});
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