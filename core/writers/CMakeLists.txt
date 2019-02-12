
add_library(gadgetron_core_writers SHARED ImageWriter.cpp ImageWriter.h BufferWriter.cpp BufferWriter.h)

target_link_libraries(gadgetron_core_writers
        gadgetron_core
        gadgetron_toolbox_mri_core)

set_target_properties(gadgetron_core_writers PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})


install(TARGETS gadgetron_core_writers DESTINATION lib COMPONENT main)
