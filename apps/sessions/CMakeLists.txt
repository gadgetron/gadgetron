include(GoogleTest)
find_package(RocksDB REQUIRED)
add_library(gadgetron_session
        SessionServer.cpp
        JSONStore.h
        RESTStorageClient.cpp
        )

target_link_libraries(gadgetron_session Boost::coroutine gadgetron_core RocksDB::RocksDB Boost::date_time )

target_include_directories(gadgetron_session PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )



add_executable(session_tests sessions_test.cpp server_test.cpp)
target_link_libraries(session_tests GTest::Main gadgetron_session gadgetron_core)

gtest_discover_tests(session_tests)