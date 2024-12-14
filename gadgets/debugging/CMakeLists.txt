include_directories(${HDF5_INCLUDE_DIRS})

set(pingvin_debugging_header_files
        PseudoReplicatorGadget.h
        RateLimitGadget.h
        WhiteNoiseInjectorGadget.h
        )

set(pingvin_debugging_src_files
        PseudoReplicatorGadget.cpp
        RateLimitGadget.cpp
        WhiteNoiseInjectorGadget.cpp
        )

add_library(pingvin_debugging SHARED
    ${pingvin_debugging_header_files}
    ${pingvin_debugging_src_files}
    )

set_target_properties(pingvin_debugging PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_debugging pingvin_core)

target_include_directories(pingvin_debugging
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(FILES
        ${pingvin_debugging_header_files}
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS pingvin_debugging
        EXPORT pingvin-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
)