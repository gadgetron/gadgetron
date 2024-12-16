add_library(pingvin_toolbox_demons SHARED demons_registration.cpp)
target_link_libraries(pingvin_toolbox_demons
    pingvin_toolbox_cpucore
    pingvin_toolbox_cpucore_math
)

set_target_properties(pingvin_toolbox_demons PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_sources(pingvin_toolbox_demons
    PUBLIC FILE_SET headers TYPE HEADERS
    FILES
        demons_registration.h
)

install(TARGETS pingvin_toolbox_demons
    EXPORT pingvin-export
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    FILE_SET headers DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH}
    INCLUDES DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH}
)
