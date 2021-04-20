
add_library(gadgetron_session
        SessionServer.cpp
        Storage.cpp
        DB.cpp
        JSONStore.h
        )

target_link_libraries(gadgetron_session Boost::coroutine gadgetron_core)

target_include_directories(gadgetron_session PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )
add_executable(gadgetron_storage StorageServer.cpp )

target_link_libraries(gadgetron_storage gadgetron_session)