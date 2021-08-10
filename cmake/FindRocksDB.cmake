
find_path(ROCKSDB_INCLUDE_DIR
        NAMES rocksdb/db.h
        )

find_library(ROCKSDB_LIBRARY NAMES rocksdb)

mark_as_advanced(ROCKSDB_FOUND ROCKSDB_INCLUDE_DIR ROCKSDB_LIBRARY)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(RocksDB REQUIRED_VARS ROCKSDB_INCLUDE_DIR ROCKSDB_LIBRARY)

if (RocksDB_FOUND)
    add_library(RocksDB::RocksDB SHARED IMPORTED)
    set_target_properties(RocksDB::RocksDB PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROCKSDB_INCLUDE_DIR}"
            IMPORTED_LOCATION ${ROCKSDB_LIBRARY})
endif()
