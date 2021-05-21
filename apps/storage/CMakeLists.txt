find_package(RocksDB REQUIRED)
add_library(gadgetron_storage
        StorageServer.cpp
        JSONStore.h
        RESTStorageClient.cpp
        )

target_link_libraries(gadgetron_storage Boost::coroutine gadgetron_core RocksDB::RocksDB Boost::date_time)

target_include_directories(gadgetron_storage PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )


add_executable(gadgetron_storage_server gadgetron_storage_server.cpp)

target_link_libraries(gadgetron_storage_server
        gadgetron_storage
        Boost::system
        Boost::filesystem
        Boost::program_options)


install(TARGETS gadgetron_storage_server DESTINATION bin COMPONENT main)

if (BUILD_TESTING)
    include(GoogleTest)
    add_executable(session_tests storage_test.cpp server_test.cpp)
    target_link_libraries(session_tests GTest::Main gadgetron_storage gadgetron_core)

    gtest_discover_tests(session_tests)
endif ()