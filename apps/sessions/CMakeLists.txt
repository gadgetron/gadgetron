include(GoogleTest)
find_package(RocksDB REQUIRED)
add_library(gadgetron_session
        SessionServer.cpp
        Storage.cpp
        DB.cpp
        JSONStore.h
        )

target_link_libraries(gadgetron_session Boost::coroutine gadgetron_core RocksDB::RocksDB Boost::date_time )

target_include_directories(gadgetron_session PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )
add_executable(gadgetron_storage StorageServer.cpp )

target_link_libraries(gadgetron_storage gadgetron_session)


add_executable(session_tests sessions_test.cpp)
target_link_libraries(session_tests GTest::Main gadgetron_session)

gtest_discover_tests(session_tests)