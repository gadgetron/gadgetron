if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_LOG__)
endif (WIN32)

add_library(gadgetron_toolbox_log SHARED log.cpp)
set_target_properties(gadgetron_toolbox_log PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_toolbox_log DESTINATION lib COMPONENT main)
install(FILES log.h log_export.h DESTINATION include COMPONENT main)
