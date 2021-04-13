
add_library(gadgetron_session
        SessionServer.cpp
        Storage.cpp
        DB.h
        )

link_libraries(gadgetron_session Boost::coroutine gadgetron_core)

add_executable(gadgetron_storage StorageServer.cpp SessionServer.cpp)

link_libraries(gadgetron_storage gadgetron_session)