#pragma once

#include <rocksdb/status.h>

namespace Gadgetron::Storage::DB {
    class DBError : public std::runtime_error {
    public:
        DBError(const rocksdb::Status &status) : std::runtime_error(status.ToString()) {}
    };
}