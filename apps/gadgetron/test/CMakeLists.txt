#include(GoogleTest)
enable_testing()

add_executable( server_tests
        socket_test.cpp ../connection/SocketStreamBuf.cpp)

target_link_libraries(server_tests
        gadgetron_core
        GTest::GTest
        GTest::Main
        gtest
        gtest_main
        )
#gtest_add_tests(TARGET server_tests AUTO)
