#pragma once


#include "Types.h"
#include <rocksdb/db.h>
#include <rocksdb/iterator.h>
#include <nlohmann/json.hpp>
#include <boost/iterator.hpp>
#include "DBError.h"
#include <range/v3/view.hpp>

namespace Gadgetron::Sessions::DB {
    using json = nlohmann::json;

    class JSONStore {
    public:
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
            auto status = database->Put(rocksdb::WriteOptions(), handle, key, str_val);
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
            using difference_type = void; //We really can't take the difference between these guys
            using pointer = const value_type *;
            using reference = value_type;

            Iterator(const Iterator &) = default;

            bool operator==(const EndPoint &endpoint) const {
                if (it->Valid()) {
                    if (endpoint.sentinel) {
                        return it->key().ToString() > endpoint.sentinel.value();
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
            friend JSONStore;
        };

        class ReverseIterator {

        public:
            using iterator_category = std::input_iterator_tag;
            using value_type = std::pair<std::string, json>;
            using difference_type = void; //We really can't take the difference between these guys
            using pointer = const value_type *;
            using reference = value_type;

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

        struct Range {
            Iterator &begin() { return front; }

            EndPoint &end() { return back; };
            Iterator front;
            EndPoint back;
        };
        struct RevRange {
            ReverseIterator &begin() { return front; }

            EndPoint &end() { return back; };
            ReverseIterator front;
            EndPoint back;
        };

       /**
         *
         * @param from Inclusive
         * @param to Inclusive
         * @return
         */
        auto range(std::string_view from, std::string_view to) {
            assert(from<to);
            auto options = rocksdb::ReadOptions();
            auto to_slice = rocksdb::Slice(to);
           auto from_slice = rocksdb::Slice(from);
           options.iterate_upper_bound = &to_slice;
            options.iterate_lower_bound = &from_slice;
            auto it = database->NewIterator(options, handle);

            return Range{Iterator(it,from),EndPoint()};
        }

        /**
        *
        * @param from Inclusive
        * @param to Inclusive
        * @return
        */
        auto reverse_range(std::string_view from, std::string_view to) {
            assert(from < to);
            auto options = rocksdb::ReadOptions();
            auto to_slice = rocksdb::Slice(to);
            options.iterate_upper_bound = &to_slice;
            auto from_slice = rocksdb::Slice(from);
            options.iterate_lower_bound = &from_slice;
            auto it = database->NewIterator(options, handle);
            return RevRange{ReverseIterator(it,to),EndPoint()};
        }

    private:

        mutable rocksdb::ColumnFamilyHandle *handle;
        mutable std::shared_ptr<rocksdb::DB> database;

    };


    template<class T>
    class ListStore {
    public:
        void push_back(std::string_view key, const T &val) {

            auto [front, back] = key_range(key);
            auto range = store.reverse_range(front,back);

            size_t front_size = front.size();

            auto get_new_index = [&]() -> uint64_t  {
                if (range.begin() == range.empty())
                    return 0;
                auto largest_key = (*range.begin()).first;
                size_t largest_index = *reinterpret_cast<size_t*>(largest_key.data()+front_size);
                return largest_index+1;

            };

            size_t index = get_new_index();
            std::string new_key(front);
            new_key.append(reinterpret_cast<const char*>(&index),sizeof(index));
            store.set(new_key,json(val));
        }

       std::vector<T> operator[](std::string_view key) {
            auto [front, back] = key_range(key);
            auto result = store.range(front, back) | ranges::view::transform([](auto&& kv){ auto& [key,value] = kv; return value.template get<T>();}) | ranges::to<std::vector>();

            return result;
        }

    private:

        std::pair<std::string, std::string> key_range(std::string_view key) const {
            std::string prefix_key(key);
            prefix_key.push_back('/');

            return {prefix_key, max_key(key)};
        }

        std::string max_key(std::string_view prefix) const {
            std::string result;
            std::uint64_t max_key_value = std::numeric_limits<std::uint64_t>::max();
            result.reserve(prefix.size() + 2 + sizeof(max_key_value));
            result.append(prefix);
            result.push_back('/');
            result.append(reinterpret_cast<const char *>(&max_key_value), sizeof(max_key_value));
            return max_key;
        }

        JSONStore store;


    };

    template<class T>
    class ValueStore {

    public:
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
};