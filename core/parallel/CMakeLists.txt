
add_library(gadgetron_core_parallel SHARED
        Branch.cpp
        Branch.h
        Merge.cpp
        Merge.h
)

target_link_libraries(gadgetron_core_parallel
        gadgetron_core
)

set_target_properties(gadgetron_core_parallel PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION}
)

install(TARGETS gadgetron_core_parallel
        DESTINATION lib COMPONENT main)

install(FILES
        Branch.h
        Merge.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
