
add_library(gadgetron_core_gadgets SHARED
        OversamplingGadget.cpp
        OversamplingGadget.h)

target_link_libraries(gadgetron_core_gadgets
        gadgetron_core
        gadgetron_toolbox_mri_core)

set_target_properties(gadgetron_core_gadgets PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_core_gadgets DESTINATION lib COMPONENT main)
