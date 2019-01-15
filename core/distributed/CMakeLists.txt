
add_library(gadgetron_core_distributed SHARED
        RemoteChannel.cpp
        RemoteChannel.h
)

target_link_libraries(gadgetron_core_distributed
        gadgetron_core
)

set_target_properties(gadgetron_core_distributed PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION}
)

install(TARGETS gadgetron_core_distributed
        DESTINATION lib COMPONENT main)

install(FILES
        RemoteChannel.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)