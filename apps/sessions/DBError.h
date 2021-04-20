#pragma once
#include <rocksdb/status.h>
class DBError : public std::runtime_error {
public:
    DBError(rocksdb::Status &status) : std::runtime_error(status.ToString()) {}
};